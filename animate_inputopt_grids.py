#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Animate grids being optimized by network

import sys, os, argparse, glob, re
import subprocess as sp
import numpy as np
from gridData import OpenDX
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

parser = argparse.ArgumentParser(description="Generate animated png of DX grids"
        " numbered by frame.")
parser.add_argument("-r", "--receptor", help="Optionally provide receptor"
        " structure(s) to include for reference in animation", 
        default='')
parser.add_argument("-l", "--ligand", help="Optionally provide ligand structure(s)"
        " to include for reference in animation", 
        default='')
parser.add_argument("-cl", "--level", default=1.0, type=float,
        help="Contour level for density maps; valid range is [-5,5], default is 1.0.")
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
args = parser.parse_args()
assert not (args.exclude_ligand and args.exclude_receptor), "Have to output at "
"least one of receptor and ligand grids"

#this is dumb but i do kind of care about mapping to reasonable colors
colormap = {}
colormap["Lig_Hydrogen"] = "hydrogen"
colormap["Lig_PolarHydrogen"] = "gray"
colormap["Lig_AliphaticCarbonXSHydrophobe"] = "lightorange"
colormap["Lig_AliphaticCarbonXSNonHydrophobe"] = "brightorange"
colormap["Lig_AromaticCarbonXSHydrophobe"] = "wheat"
colormap["Lig_AromaticCarbonXSNonHydrophobe"] = "brown"
colormap["Lig_Nitrogen"] = "deepblue"
colormap["Lig_NitrogenXSDonor"] = "br3"
colormap["Lig_NitrogenXSDonorAcceptor"] = "br4"
colormap["Lig_NitrogenXSAcceptor"] = "br5"
colormap["Lig_Oxygen"] = "firebrick"
colormap["Lig_OxygenXSDonor"] = "br7"
colormap["Lig_OxygenXSDonorAcceptor"] = "br8"
colormap["Lig_OxygenXSAcceptor"] = "br9"
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
colormap["Rec_AliphaticCarbonXSHydrophobe"] = "lightorange"
colormap["Rec_AliphaticCarbonXSNonHydrophobe"] = "brightorange"
colormap["Rec_AromaticCarbonXSHydrophobe"] = "wheat"
colormap["Rec_AromaticCarbonXSNonHydrophobe"] = "brown"
colormap["Rec_Nitrogen"] = "deepblue"
colormap["Rec_NitrogenXSDonor"] = "br3"
colormap["Rec_NitrogenXSDonorAcceptor"] = "br4"
colormap["Rec_NitrogenXSAcceptor"] = "br5"
colormap["Rec_Oxygen"] = "firebrick"
colormap["Rec_OxygenXSDonor"] = "br7"
colormap["Rec_OxygenXSDonorAcceptor"] = "br8"
colormap["Rec_OxygenXSAcceptor"] = "br9"
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

#figure out the grids; make a dict with a key for each channel mapped to a
#list of the frames in order
frames = {}
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
for g in grids:
    match = re.match(pattern,g)
    if match:
        gridclass = match.group(3)
        if args.exclude_receptor and gridclass == "Rec":
            continue
        elif args.exclude_ligand and gridclass == "Lig":
            continue
        channel = '%s_%s' %(gridclass, match.group(4))
        if channel not in frames:
            frames[channel] = []
        framenum = int(match.group(2).lstrip('iter'))
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

print "Generating images\n"
with open('animate.py', 'w') as f:
    f.write('from pymol import cmd\n\n')
    for i in range(maxframes+1):
        if args.receptor:
            f.write('cmd.load("%s", "rec")\n' %args.receptor)
            f.write('cmd.hide("lines", "rec")\n')
            f.write('cmd.cartoon("automatic", "rec")\n')
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
                    raise FileNotFoundError('%s not found\n' %fname)
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
