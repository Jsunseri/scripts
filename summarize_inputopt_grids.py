#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Summarize grids produced by input optimization. Currently supports producing
# an animation of the optimization process, summarizing average grid changes
# over time (mostly to assess convergence), showing how much grids changed over
# time at a per-voxel granularity (suitable e.g. for animating the figures)

import sys, os, argparse, glob, re, math, cPickle
import subprocess as sp
from multiprocessing import Pool, RawArray
from functools import partial
import numpy as np
from gridData import OpenDX
import matplotlib.pyplot as plt
import matplotlib as mpl
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()
mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)

def debug_signal_handler(signal, frame):
    import pdb
    pdb.set_trace()
import signal
signal.signal(signal.SIGINT, debug_signal_handler)

def get_plot_columns(channel_list, channel_groups):
    '''
    Given the channels being used and the way we want to group the columns,
    figure out the number of rows and columns as well as the mapping from
    channel and plot id in the final plot grid
    '''
    plot_ids = {}
    plot_groups = [[] for _ in xrange(len(channel_groups))]
    for i,channels in enumerate(channel_groups):
        plot_groups[i] = [el for el in channels if el in channel_list]
    plot_groups = [olist for olist in plot_groups if olist]
    grid_length = max([len(l) for l in plot_groups])
    grid_width = len(plot_groups)
    for i,column in enumerate(plot_groups):
        for j,channel in enumerate(column):
            plot_ids[channel] = i + j * grid_width
    return (grid_length, grid_width, plot_ids)

plot_options = ['movie', 'convergence', 'per_voxel', 'images']

class DefaultListAction(argparse.Action):
    '''
    Support providing list of arguments with a default list that is replaced
    with the user selection rather than appended to
    '''
    options = plot_options
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            for value in values:
                if value not in self.options:
                    message = ("invalid choice: {0!r} (choose from"
                            " {1})".format(value, ', '.join([repr(action) for
                                action in self.options])))
                    raise argparse.ArgumentError(self, message)
            setattr(namespace, self.dest, values)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate summary visualizations"
            "after performing input optimization, including a movie of grids optimizing,"
            "per-iteration changes in average density, and overall changes per-voxel.")
    parser.add_argument("-r", "--receptor", help="Optionally provide receptor"
            " structure(s) to include for reference in animation", 
            default='')
    parser.add_argument("-l", "--ligand", help="Optionally provide ligand structure(s)"
            " to include for reference in animation", 
            default='')
    parser.add_argument("-cl", "--level", default=1.0, type=float,
            help="Contour level for density maps; valid range is [-5,5], default is 1.0.")
    parser.add_argument("-sn", "--show_negative", default=False,
            action="store_true", help="Show negative regions of density maps as"
            "isomesh of the same color as the positive density")
    parser.add_argument("-m", "--maps", help="Optionally provide density maps to"
            " use; if not provided we'll glob for them in the current directory", 
            nargs='*', default=[])
    parser.add_argument("-t", "--threshold", default=0.0, type=float, help="Specify"
            " threshold for delta between start and end densities, below which a"
            " density will not be included")
    parser.add_argument("-er", "--exclude_receptor", default=False,
            action="store_true", help="Exclude receptor density")
    parser.add_argument("-el", "--exclude_ligand", default=False,
            action="store_true", help="Exclude ligand density")
    parser.add_argument("--sigma", default=False, action="store_true",
            help="Display sigma value used for isosurface on final animation")
    parser.add_argument("-v", "--view", default='', help='Optionally'
            ' provide a view for pymol for stricter control over the final scene,'
            ' must be formatted as a string of eighteen comma-separated floats')
    parser.add_argument("-k", "--keep_start", default=False, action='store_true',
            help='Keep whatever density is present at iteration 0 for reference')
    parser.add_argument("--figs", nargs="*", action=DefaultListAction,
            default=plot_options, metavar='ACTION')
    parser.add_argument("--pickle", default='', help='Provide pickle of'
            ' frames dict mapping desired channels to frames, per-channel deltas'
            ' per timestep')
    parser.add_argument("-O", "--overwrite", default=False, action='store_true',
            help='Overwrite existing image files in this directory')
    args = parser.parse_args()
    assert not (args.exclude_ligand and args.exclude_receptor), "Have to output at "
    "least one of receptor and ligand grids"
    
    #this is dumb but i do kind of care about mapping to reasonable colors
    colormap = {}
    colormap["Lig_Hydrogen"] = "hydrogen"
    colormap["Lig_PolarHydrogen"] = "gray"
    colormap["Lig_AliphaticCarbonXSHydrophobe"] = "orange"
    colormap["Lig_AliphaticCarbonXSNonHydrophobe"] = "o2"
    colormap["Lig_AromaticCarbonXSHydrophobe"] = "brown"
    colormap["Lig_AromaticCarbonXSNonHydrophobe"] = "wheat"
    colormap["Lig_Nitrogen"] = "lightblue"
    colormap["Lig_NitrogenXSDonor"] = "deepblue"
    colormap["Lig_NitrogenXSDonorAcceptor"] = "marine"
    colormap["Lig_NitrogenXSAcceptor"] = "purpleblue"
    colormap["Lig_Oxygen"] = "salmon"
    colormap["Lig_OxygenXSDonor"] = "raspberry"
    colormap["Lig_OxygenXSDonorAcceptor"] = "deepsalmon"
    colormap["Lig_OxygenXSAcceptor"] = "tv_red"
    colormap["Lig_Sulfur"] = "sulfur"
    colormap["Lig_SulfurAcceptor"] = "olive"
    colormap["Lig_Phosphorus"] = "hotpink"
    colormap["Lig_Fluorine"] = "limegreen"
    colormap["Lig_Chlorine"] = "lime"
    colormap["Lig_Bromine"] = "limon"
    colormap["Lig_Iodine"] = "palecyan"
    colormap["Lig_Magnesium"] = "silver"
    colormap["Lig_Manganese"] = "dirtyviolet"
    colormap["Lig_Zinc"] = "zinc"
    colormap["Lig_Calcium"] = "calcium"
    colormap["Lig_Iron"] = "iron"
    colormap["Lig_GenericMetal"] = "lightblue"
    colormap["Lig_Boron"] = "blue"
    colormap["Rec_Hydrogen"] = "hydrogen"
    colormap["Rec_PolarHydrogen"] = "gray"
    colormap["Rec_AliphaticCarbonXSHydrophobe"] = "orange"
    colormap["Rec_AliphaticCarbonXSNonHydrophobe"] = "o2"
    colormap["Rec_AromaticCarbonXSHydrophobe"] = "brown"
    colormap["Rec_AromaticCarbonXSNonHydrophobe"] = "wheat"
    colormap["Rec_Nitrogen"] = "lightblue"
    colormap["Rec_NitrogenXSDonor"] = "deepblue"
    colormap["Rec_NitrogenXSDonorAcceptor"] = "marine"
    colormap["Rec_NitrogenXSAcceptor"] = "purpleblue"
    colormap["Rec_Oxygen"] = "salmon"
    colormap["Rec_OxygenXSDonor"] = "raspberry"
    colormap["Rec_OxygenXSDonorAcceptor"] = "deepsalmon"
    colormap["Rec_OxygenXSAcceptor"] = "tv_red"
    colormap["Rec_Sulfur"] = "sulfur"
    colormap["Rec_SulfurAcceptor"] = "olive"
    colormap["Rec_Phosphorus"] = "hotpink"
    colormap["Rec_Fluorine"] = "limegreen"
    colormap["Rec_Chlorine"] = "lime"
    colormap["Rec_Bromine"] = "limon"
    colormap["Rec_Iodine"] = "palecyan"
    colormap["Rec_Magnesium"] = "silver"
    colormap["Rec_Manganese"] = "dirtyviolet"
    colormap["Rec_Zinc"] = "zinc"
    colormap["Rec_Calcium"] = "calcium"
    colormap["Rec_Iron"] = "iron"
    colormap["Rec_GenericMetal"] = "lightblue"
    colormap["Rec_Boron"] = "blue"
    
    column_groups = []
    column_groups.append(['Lig_AliphaticCarbonXSHydrophobe',
            'Lig_AliphaticCarbonXSNonHydrophobe', 'Lig_AromaticCarbonXSHydrophobe', 
            'Lig_AromaticCarbonXSNonHydrophobe'])
    column_groups.append(['Lig_Nitrogen', 'Lig_NitrogenXSDonor', 'Lig_NitrogenXSDonorAcceptor',
            'Lig_NitrogenXSAcceptor'])
    column_groups.append(['Lig_Oxygen', 'Lig_OxygenXSDonor', 'Lig_OxygenXSDonorAcceptor',
        'Lig_OxygenXSAcceptor'])
    column_groups.append(['Lig_Fluorine', 'Lig_Chlorine', 'Lig_Bromine', 'Lig_Iodine'])
    column_groups.append(['Lig_Sulfur', 'Lig_SulfurAcceptor', 'Lig_Phosphorus',
        'Lig_Boron'])
    column_groups.append(['Lig_Magnesium', 'Lig_Manganese', 'Lig_Zinc', 'Lig_Calcium'])
    column_groups.append(['Lig_Iron', 'Lig_GenericMetal'])
    
    frames = {}
    cdiffs = []
    ymax = 0
    ymin = 0
    #if pickle was provided, sanity check it
    if args.pickle:
        if not os.path.isfile(args.pickle):
            raise FileNotFoundError('%s not found\n' %g)
        frames,cdiffs,ymin,ymax = cPickle.load(open(args.pickle))
        if frames:
            try:
                maxframes = max([len(iters) for iters in frames.values()])-1
            except (AttributeError, TypeError):
                raise AssertionError('frames should be a dict mapping channels to frames')
            for channel,filelist in frames.iteritems():
                try:
                    colormap[channel]
                except KeyError:
                    raise AssertionError('frames have unknown channel %s'
                            %channel)
                for (f,num) in filelist:
                    if not os.path.isfile(f):
                        raise FileNotFoundError('%s not found\n' %f)
                try:
                    frames[channel].sort(key=lambda x:x[1])
                except:
                    raise AssertionError('frames should have lists of file:framenum tuples')
    if not frames:
        #figure out the grids; make a dict with a key for each channel mapped to a
        #list of the frames in order
        maxframes = 0
        #if list of grids was provided
        #if not, we glob for dx files 
        grids = args.maps
        if not grids:
            grids = glob.glob("*.dx")
        else:
            for g in grids:
                if not os.path.isfile(g):
                    raise FileNotFoundError('%s not found\n' %g)
        
        #attempt to parse the file names according to what the input opt script will
        #name them, i.e. [recname]_[ligname]_[iternum]_[channel].dx
        pattern = r'([A-Za-z\d_]+)_(iter\d+)_(Rec|Lig)_([A-Za-z]+).dx$'
        alt_pattern = r'((iter\d+)_[A-Za-z\d_]+)_(Rec|Lig)_([A-Za-z]+).dx$'
        for g in grids:
            base = os.path.basename(g)
            match = re.match(pattern,base)
            altmatch = re.match(alt_pattern,base)
            if match or altmatch:
                thismatch = match if match else altmatch
                gridclass = thismatch.group(3)
                if args.exclude_receptor and gridclass == "Rec":
                    continue
                elif args.exclude_ligand and gridclass == "Lig":
                    continue
                channel = '%s_%s' %(gridclass, thismatch.group(4))
                if channel not in frames:
                    frames[channel] = []
                if match:
                    framenum = int(match.group(2).lstrip('iter'))
                elif altmatch:
                    framenum = int(''.join([char for char in altmatch.group(1) if
                        char.isdigit()]))
                if framenum > maxframes: maxframes = framenum
                frames[channel].append((g,framenum))
            else:
                print "%s could not be matched to pattern and will be excluded. Continuing...\n" %g
        
        if not frames:
            print "No grid names match [recname]_[ligname]_[iternum]_[channel].dx pattern, so I can't parse them. Exiting...\n"
            sys.exit()
        
        for channel in frames:
            frames[channel].sort(key=lambda x:x[1])
        
        print "Got files...\n"
        #if we have a delta threshold, assess each channel to see if it meets the
        #threshold and pop any that don't
        #some question of how to apply threshold - sum up all the voxels? calculate
        #mean intensity? for now, looking at maximum absolute change per voxel
        #TODO: parallelize this or otherwise speed it up
        frame0 = {}
        origin = np.array([0, 0, 0])
        if args.threshold:
            print "Applying threshold...\n"
            for channel in frames.keys():
                start = OpenDX.field(0)
                start.read(frames[channel][0][0])
                origin = start.components['positions'].origin
                sdat = start.components['data'].array
                if args.keep_start:
                    frame0[channel] = frames[channel][0][0]
                sufficient = False
                framenums = range(len(frames[channel]))
                framenums = [framenums[-1-i] for i in range(len(framenums))]
                #for now
                framenums = [len(frames[channel])-1]
                for idx in framenums:
                    print "%s,%d\n" %(channel,idx)
                    nextframe = OpenDX.field(0)
                    nextframe.read(frames[channel][idx][0])
                    ndat = nextframe.components['data'].array
                    if np.max(np.abs(sdat - ndat)) > args.threshold:
                        sufficient = True
                        break
                if not sufficient:
                    frames.pop(channel, None)
        
        if not frames:
            print "No channels changed enough to pass the threshold. Exiting.\n"
            sys.exit()

    if 'images' in args.figs:
        print "Generating images\n"
        with open('make_frames.py', 'w') as f:
            f.write('from pymol import cmd\n\n')
            #take the length of the list?
            firstkey = frames.keys()[0]
            numframes = len(frames[firstkey])
            for i in range(numframes):
                outname = 'frame%d.png' %i
                if not args.overwrite and os.path.isfile(outname):
                    continue
                if args.receptor:
                    f.write('cmd.load("%s", "rec")\n' %args.receptor)
                    f.write('cmd.hide("lines", "rec")\n')
                    f.write('cmd.show("cartoon", "rec")\n')
                if args.ligand:
                    f.write('cmd.load("%s", "lig")\n' %args.ligand)
                    f.write('cmd.show("sticks", "lig")\n')
                f.write('cmd.bg_color("white")\n')
                f.write('util.cbaw()\n')
                f.write('cmd.set("ray_opaque_background", "off")\n')
                f.write('cmd.set("ray_shadows", "off")\n')
                f.write('cmd.set("ray_volume", "on")\n')
                f.write('cmd.set("antialias", 2)\n')
                f.write('cmd.set_color("o2",[1,0.75,0])\n')
                if args.show_negative:
                    f.write('cmd.set("mesh_negative_visible")\n')
                for channel in frames:
                    if len(frames[channel]) >= (i+1):
                        if channel.startswith('Rec'):
                            level = 1.0
                        else:
                            level = args.level
                        fname = frames[channel][i][0]
                        if not os.path.isfile(fname):
                            raise OSError('%s not found\n' %fname)
                        f.write('cmd.load("%s")\n' %fname)
                        objname = "%s_%d" %(channel,i)
                        if args.show_negative:
                            f.write('cmd.isomesh("%s", "%s", level=%f)\n'
                                    %(objname+'_mesh',
                                        os.path.splitext(fname)[0], level))
                            f.write('cmd.color("%s", "%s")\n'
                                    %(colormap[channel],objname+'_mesh'))
                            f.write('cmd.set("mesh_negative_color", "%s", "%s")\n'
                                    %(colormap[channel],objname+'_mesh'))
                        f.write('cmd.isosurface("%s", "%s", level=%f)\n'
                                %(objname,os.path.splitext(fname)[0],level))
                        f.write('cmd.color("%s", "%s")\n' %(colormap[channel],objname))
                if not args.view:
                    f.write('cmd.origin(position=[%f, %f, %f])\n' %(origin[0], origin[1],
                        origin[2]))
                    f.write('cmd.center("origin")\n')
                    f.write('cmd.orient("lig", state=1)\n')
                    f.write('cmd.zoom("center", 20, complete=1)\n')
                else:
                    f.write('cmd.set_view("%s")\n' %args.view)
                f.write('cmd.png("frame%d_level%f.png", width=1080, height=1080, dpi=300, ray=1)\n' %(i, args.level))
                f.write('cmd.save("frame%d_level%f.pse")\n' %(i,args.level))
                f.write('cmd.reinitialize()\n')
        cmd = 'pymol -c make_frames.py'
        p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
        out,err = p.communicate()
        if err:
            print err

    #we processed them the first time to get the frames in order; now we add
    #padding (empty grids) where necessary to deal with any gaps
    final_frames = {}
    for channel in frames:
        final_frames[channel] = []
        for element in frames[channel]:
            true_num = element[1]
            while (true_num > len(final_frames[channel])):
                next_num = len(final_frames[channel])
                final_frames[channel].append(('', next_num))
            final_frames[channel].append(element)
    
    frames = final_frames
    
    if 'movie' in args.figs:
        print "Generating movie images\n"
        with open('animate.py', 'w') as f:
            f.write('from pymol import cmd\n\n')
            for i in range(maxframes+1):
                outname = 'frame%d.png' %i
                if not args.overwrite and os.path.isfile(outname):
                    continue
                if args.receptor:
                    f.write('cmd.load("%s", "rec")\n' %args.receptor)
                    f.write('cmd.hide("lines", "rec")\n')
                    f.write('cmd.show("cartoon", "rec")\n')
                if args.ligand:
                    f.write('cmd.load("%s", "lig")\n' %args.ligand)
                    f.write('cmd.show("sticks", "lig")\n')
                f.write('cmd.bg_color("white")\n')
                f.write('util.cbaw()\n')
                f.write('cmd.set("ray_opaque_background", "off")\n')
                f.write('cmd.set("ray_shadows", "off")\n')
                f.write('cmd.set("ray_volume", "on")\n')
                f.write('cmd.set("antialias", 2)\n')
                if args.keep_start:
                    for channel in frame0:
                        fname = frame0[channel]
                        f.write('cmd.load("%s")\n' %fname)
                        objname = "%s_%s" %(channel,'start')
                        f.write('cmd.isosurface("%s", "%s", level=%f)\n'
                                %(objname,os.path.splitext(fname)[0],1))
                        f.write('cmd.color("%s", "%s")\n' %(colormap[channel],objname))
                for channel in frames:
                    assert frames[channel][i][1] == i, "framenum %s at index %d" %(frames[name][channel][i][1],i)
                    if len(frames[channel]) >= (i+1):
                        fname = frames[channel][i][0]
                        if not os.path.isfile(fname):
                            raise OSError('%s not found\n' %fname)
                        f.write('cmd.load("%s")\n' %fname)
                        objname = "%s_%d" %(channel,i)
                        f.write('cmd.isosurface("%s", "%s", level=%f)\n'
                                %(objname,os.path.splitext(fname)[0],args.level))
                        f.write('cmd.color("%s", "%s")\n' %(colormap[channel],objname))
                if not args.view:
                    f.write('cmd.origin(position=[%f, %f, %f])\n' %(origin[0], origin[1],
                        origin[2]))
                    f.write('cmd.center("origin")\n')
                    f.write('cmd.orient("lig", state=1)\n')
                    f.write('cmd.zoom("center", 20, complete=1)\n')
                else:
                    f.write('cmd.set_view("%s")\n' %args.view)
                f.write('cmd.png("frame%d.png", width=1080, height=1080, dpi=300, ray=1)\n' %(i))
                f.write('cmd.reinitialize()\n')
        cmd = 'pymol -c animate.py'
        p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
        out,err = p.communicate()
        if err:
            print err
        
        print "Making movie\n"
        if args.sigma:
            cmd = 'ffmpeg -y -framerate 5 -i %s%%d.png -plays 0 -vf '
            'drawtext="fontfile=/usr/share/fonts/truetype/cmu/cmuntt.ttf: '
            'text=\'σ=%.2f\':fontsize=24" %s.apng' %('frame', args.level, 'frame')
        else:
            cmd = 'ffmpeg -y -framerate 5 -i %s%%d.png -plays 0 %s.mp4' %('frame',
                    'frame')
        p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
        out,err = p.communicate()
        if err:
            print err
    
    def pairwise_timestep_diffs(channel_frames):
        '''
        Returns list of mean diffs per timestep for a given channel
        '''
        maxframes = len(channel_frames)
        cdiffs = []
        t1 = OpenDX.field(0)
        t1.read(channel_frames[0][0])
        t1_pixels = t1.components['data'].array
        for t in xrange(1, maxframes):
            t2 = OpenDX.field(0)
            t2.read(channel_frames[t][0])
            t2_pixels = t2.components['data'].array
            mean_diff = np.mean(t2_pixels - t1_pixels)
            cdiffs.append(mean_diff)
            t1_pixels = t2_pixels
        return cdiffs
    
    if 'convergence' in args.figs:
        #if convergence, make per-channel and overall plots of average T_i - T_i-1 per-voxel over time
        print "Summarizing per-voxel change over time"
        if cdiffs:
            #got them from the pickle, let's sanity check
            assert len(cdiffs) == len(frames), "Different numbers of channels in cdiffs and frames"
            nsteps = len(frames.values()[0])
            for i,tsteps in enumerate(cdiffs):
                this_tsteps = len(tsteps)
                assert this_tsteps == nsteps, "Incorrect number of timesteps at %d, %d vs %d\n" %(i, this_tsteps, nsteps)
        else:
            pool = Pool()
            cdiffs = pool.map(pairwise_timestep_diffs, frames.values())
    
        #average over all channels
        fig,ax = plt.subplots(figsize=(16,16), num=1)
        #per-channel
        cfig,cax = plt.subplots(figsize=(16,16), num=2)
        for i,channel in enumerate(frames.iterkeys()):
            cax.plot(range(maxframes), cdiffs[i],
                    color=pymol.querying.get_color_tuple(colormap[channel]),
                    label=channel.replace('_', ' '), lw=6)
    
        #per-channel fig
        cax.set_xlabel('Iteration')
        cax.set_ylabel('Average per-pixel change')
        cax.legend(loc='best', frameon=True)
        cfig.savefig('per_channel_change.pdf', bbox_inches='tight')
    
        #mean over all channels
        diffs = np.mean(np.array(cdiffs), axis=0)
        ax.plot(range(maxframes), diffs, lw=6)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Average per-pixel change')
        fig.savefig('overall_change.pdf', bbox_inches='tight')
    
    def cumulative_voxeldiff_plot(frames, plot_ids, grid_length, grid_width,
            start_data, ymin, ymax, t=0, force=False):
        '''
        for timestep t, compute per-voxel change since the start of optimization
        and plot it
        '''
        outname = "voxels_iter%d.png" %t
        if not args.overwrite and not force and os.path.isfile(outname):
            return (ymin, ymax)
        print "Plotting timestep %d" %t
        fig,ax = plt.subplots(figsize=(16,16))
        channel_ax = []
        for channel in frames:
            plot_num = plot_ids[channel]
            tnow = OpenDX.field(0)
            tnow.read(frames[channel][t][0])
            tnow_pixels = tnow.components['data'].array
            subax = plt.subplot2grid((grid_length, grid_width),
                (plot_num/grid_width, plot_num % grid_width), fig=fig)
            subax.plot(range(tnow_pixels.size), tnow_pixels - start_data[channel],
                    color=pymol.querying.get_color_tuple(colormap[channel]),
                    ls=None)
            subax.set_title(channel.replace('_', ' '))
            if plot_num % grid_width == 0:
                subax.set_ylabel('Density')
            if plot_num / grid_width == grid_length-1:
                subax.set_xlabel('Voxel')
            for tick in subax.get_xticklabels():
                tick.set_rotation(45)
            channel_ax.append(subax)
            ylims = subax.get_ylim()
            if ylims[0] < ymin:
                ymin = ylims[0]
            if ylims[1] > ymax:
                ymax = ylims[1]
        if ymax:
            for i,_ in enumerate(frames):
                channel_ax[i].set_ylim(ymin, ymax)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        fig.savefig(outname, bbox_inches="tight")
        plt.close(fig)
        return (ymin,ymax)
    
    if 'per_voxel' in args.figs:
        #if per_voxel, for each timestep generate a grid of plots (one for each
        #channel) showing actual change per voxel compared with the start
        #and let's animate it just cuz
        print "Explicitly plotting per-voxel changes, per-channel, over time"
        grid_length,grid_width,plot_ids = get_plot_columns(frames.keys(),
                column_groups)
        start_data = {}
        for channel in frames:
            start = OpenDX.field(0)
            start.read(frames[channel][0][0])
            start_data[channel] = start.components['data'].array
    
        if not ymax:
            #if not ymax, we didn't get it from a pickle; let's try to set it
            #reasonably at the start by doing the last frame first before chunking
            ymin,ymax = cumulative_voxeldiff_plot(frames, plot_ids, grid_length, grid_width,
                    start_data, ymin, ymax, maxframes, True)
    
        #originally went backwards to set ymin/ymax for the last timestep
        pool = Pool(4)
        func = partial(cumulative_voxeldiff_plot, frames, plot_ids,
                grid_length, grid_width, start_data, ymin, ymax)
        framenums = range(maxframes-1, 0, -1)
        ylims = pool.map(func, framenums)
        real_ymin = min(ylims, lambda x: x[0])
        real_ymax = max(ylims, lambda x: x[1])
        if real_ymin < ymin:
            print "Inconcistent ymin over frames\n"
        if real_ymax > ymax:
            print "Inconcistent ymax over frames\n"
    
        cmd = 'ffmpeg -y -framerate 5 -i %s%%d.png -plays 0 %s.mp4' %('voxels_iter',
                'voxels')
        p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
        out,err = p.communicate()
        if err:
            print err
    
    cPickle.dump((frames,cdiffs,ymin,ymax), open('summary.cpickle', 'w'), -1)
