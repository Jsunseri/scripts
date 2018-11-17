#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Animate grids being optimized by network

import os, argparse, glob, re
import subprocess as sp
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

parser = argparse.ArgumentParser(description="Generate animated png of DX grids"
        " numbered by frame.")
parser.add_argument("-r", "--receptor", help="Optionally provide receptor"
        " structure(s) to include for reference in animation, frames and"
        " structure should have a shared name prefix that can be extracted by"
        " splitting on underscores", nargs='*', default=[])
parser.add_argument("-l", "--ligand", help="Optionally provide ligand structure(s)"
        " to include for reference in animation, frames and structure should"
        " have a shared name prefix that can be extracted by splitting on"
        " underscores", nargs='*', default=[])
parser.add_argument("-cl", "--level", default=1.0, type=float,
        help="Contour level for density maps; valid range is [-5,5], default is 1.0.")
parser.add_argument("-er", "--exclude_receptor", default=False,
        action="store_true")
parser.add_argument("-el", "--exclude_ligand", default=False,
        action="store_true")
args = parser.parse_args()
assert not (args.exclude_ligand and args.exclude_receptor), "Have to output at "
"least one of receptor and ligand grids"

#this is dumb but i do kind of care about mapping to reasonable colors
colormap = {}
colormap["lig_Hydrogen"] = "hydrogen"
colormap["lig_PolarHydrogen"] = "gray"
colormap["lig_AliphaticCarbonXSHydrophobe"] = "lightorange"
colormap["lig_AliphaticCarbonXSNonHydrophobe"] = "brightorange"
colormap["lig_AromaticCarbonXSHydrophobe"] = "wheat"
colormap["lig_AromaticCarbonXSNonHydrophobe"] = "brown"
colormap["lig_Nitrogen"] = "deepblue"
colormap["lig_NitrogenXSDonor"] = "br3"
colormap["lig_NitrogenXSDonorAcceptor"] = "br4"
colormap["lig_NitrogenXSAcceptor"] = "br5"
colormap["lig_Oxygen"] = "firebrick"
colormap["lig_OxygenXSDonor"] = "br7"
colormap["lig_OxygenXSDonorAcceptor"] = "br8"
colormap["lig_OxygenXSAcceptor"] = "br9"
colormap["lig_Sulfur"] = "sulfur"
colormap["lig_SulfurAcceptor"] = "olive"
colormap["lig_Phosphorus"] = "hotpink"
colormap["lig_Fluorine"] = "limegreen"
colormap["lig_Chlorine"] = "lime"
colormap["lig_Bromine"] = "limon"
colormap["lig_Iodine"] = "palecyan"
colormap["lig_Magnesium"] = "silver"
colormap["lig_Manganese"] = "dirtyviolet"
colormap["lig_Zinc"] = "zinc"
colormap["lig_Calcium"] = "calcium"
colormap["lig_Iron"] = "iron"
colormap["lig_GenericMetal"] = "lightblue"
colormap["lig_Boron"] = "blue"
colormap["rec_Hydrogen"] = "hydrogen"
colormap["rec_PolarHydrogen"] = "gray"
colormap["rec_AliphaticCarbonXSHydrophobe"] = "lightorange"
colormap["rec_AliphaticCarbonXSNonHydrophobe"] = "brightorange"
colormap["rec_AromaticCarbonXSHydrophobe"] = "wheat"
colormap["rec_AromaticCarbonXSNonHydrophobe"] = "brown"
colormap["rec_Nitrogen"] = "deepblue"
colormap["rec_NitrogenXSDonor"] = "br3"
colormap["rec_NitrogenXSDonorAcceptor"] = "br4"
colormap["rec_NitrogenXSAcceptor"] = "br5"
colormap["rec_Oxygen"] = "firebrick"
colormap["rec_OxygenXSDonor"] = "br7"
colormap["rec_OxygenXSDonorAcceptor"] = "br8"
colormap["rec_OxygenXSAcceptor"] = "br9"
colormap["rec_Sulfur"] = "sulfur"
colormap["rec_SulfurAcceptor"] = "olive"
colormap["rec_Phosphorus"] = "hotpink"
colormap["rec_Fluorine"] = "limegreen"
colormap["rec_Chlorine"] = "lime"
colormap["rec_Bromine"] = "limon"
colormap["rec_Iodine"] = "palecyan"
colormap["rec_Magnesium"] = "silver"
colormap["rec_Manganese"] = "dirtyviolet"
colormap["rec_Zinc"] = "zinc"
colormap["rec_Calcium"] = "calcium"
colormap["rec_Iron"] = "iron"
colormap["rec_GenericMetal"] = "lightblue"
colormap["rec_Boron"] = "blue"

#figure out the grids; make a dict with a key for each channel mapped to a
#list of the frames in order
frames = {}
maxframes = 0
grids = glob.glob("*.dx")
pattern = r'([\w\.]+?)gninatypes_(rec|lig)_([A-Za-z]+)(\d+)\.dx$'
for g in grids:
    name = g.split('_')[0]
    if name not in frames:
        frames[name] = {}
    match = re.match(pattern,g)
    if match:
        gridclass = match.group(2)
        if args.exclude_receptor and gridclass == "rec":
            continue
        elif args.exclude_ligand and gridclass == "lig":
            continue
        channel = '%s_%s' %(gridclass, match.group(3))
        if channel not in frames[name]:
            frames[name][channel] = []
        framenum = int(match.group(4))
        if framenum > maxframes: maxframes = framenum
        frames[name][channel].append((g,framenum))

for name in frames:
    for channel in frames[name]:
        frames[name][channel].sort(key=lambda x:x[1])

#we processed them the first time to get the frames in order; now we add
#padding where necessary to deal with gaps
#yuck
final_frames = {}
for name in frames:
    final_frames[name] = {}
    for channel in frames[name]:
        final_frames[name][channel] = []
        for num in range(len(frames[name][channel])):
            true_num = frames[name][channel][num][1]
            while (true_num > len(final_frames[name][channel])):
                final_frames[name][channel].append(('', ''))
            final_frames[name][channel].append(frames[name][channel][num])

frames = final_frames
#if we have rec/lig structures, we attempt to match them by "name"...could make this more flexible
recs = {}
for rec in args.receptor:
    recs[os.path.basename(rec).split('_')[0]] = rec
ligs = {}
for lig in args.ligand:
    ligs[os.path.basename(lig).split('_')[0]] = lig

for name in frames:
    if name in recs:
        pymol.cmd.load(recs[name], "rec")
        pymol.cmd.hide("lines", "rec")
        pymol.cmd.cartoon("automatic", "rec")
    if name in ligs:
        pymol.cmd.load(ligs[name], "lig")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_shadows", "off")
    pymol.cmd.set("antialias", 2)
    for i in range(maxframes+1):
        for channel in frames[name]:
            if len(frames[name][channel]) >= (i+1) and frames[name][channel][i][0]:
                fname = frames[name][channel][i][0]
                assert frames[name][channel][i][1] == i, "framenum %s at index %d" %(frames[name][channel][i][1],i)
                pymol.cmd.load(fname)
                objname = "%s_%d" %(channel,i)
                pymol.cmd.isosurface(objname,os.path.splitext(fname)[0],args.level)
                pymol.cmd.color(colormap[channel],objname)
        pymol.cmd.png("%s_%d.png" %(name, i), width=1080, height=1080, dpi=100, ray=1)
    pymol.cmd.reinitialize()

    cmd = 'ffmpeg -y -framerate 2 -i %s_%%d.png -plays 0 -vf '
    'drawtext="fontfile=/usr/share/fonts/truetype/cmu/cmuntt.ttf: '
    'text=\'σ=%.2f\':fontsize=24" %s.apng' %(name, args.level, name)
    p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
    out,err = p.communicate()
    if err:
            print err
