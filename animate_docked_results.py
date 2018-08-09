#!/usr/bin/python

# Creates animated GIF cycling through docked poses

import sys, os, math, random
from rdkit.Chem import AllChem
import subprocess as sp
import argparse

def count_mols(file):
    count = 0
    with open(file) as f:
        for line in f:
            if line.startswith('$$$$'):
                count += 1
    return count

def make_directory(base):
	""" 
	Make directory with error checking that 
	avoids overwriting an existing directory.
	"""
		
	i = 0
	while 1:
	    try:
	        if i == 0:
	            dirname = base
	        else:
	            dirname = base + '_' + str(i)
	        os.mkdir(dirname)
	        break
	    except OSError:
	        if not os.path.isdir(dirname):
	            raise
	        i += 1
	        pass 
	return dirname

parser = argparse.ArgumentParser(description="Generates gif of docked \
poses of input small molecule, optionally in a binding site \
(prealignment assumed).")
parser.add_argument('-l', '--ligand', help='Small molecule whose provided poses \
will be sampled')
parser.add_argument('-r', '--receptor', default='', help='Optionally \
include receptor for animation')
parser.add_argument('-o', '--out', default='out', help='Specify outprefix, \
default is "out"')

args = parser.parse_args()
ext = os.path.splitext(args.ligand)[1]
if ext != ".sdf":
    print "Only SDFs supported for ligand"
    sys.exit()

num_mols = count_mols(args.ligand)
dirname = make_directory(args.out)
os.chdir(dirname)
with open('friends.py','w') as file:
    file.write('from pymol import cmd\n\n')
    file.write('cmd.load("%s")\n' %args.ligand)
    file.write('cmd.bg_color("white")\n')
    file.write('cmd.show_as("sticks")\n')
    file.write('util.cbaw()\n')
    file.write('cmd.set("ray_opaque_background", "off")\n')
    file.write('cmd.set("ray_shadows", "off")\n')
    file.write('cmd.set("antialias", 2)\n')
    file.write('cmd.orient("all", state=1)\n')
    #file.write('cmd.set("ray_trace_frames","1")\n')
    if args.receptor:
        file.write('cmd.load("%s")\n' %args.receptor)
        recname = os.path.splitext(os.path.basename(args.receptor))[0]
        file.write('cmd.show("surface", "%s")\n' %recname)
        file.write('cmd.set("transparency", "0.2")\n')
        file.write('cmd.set("stick_radius", "0.15", "%s")\n' %recname)
    for i in range(1,num_mols+1):
        file.write('cmd.set("state",%d)\n' %i)
        file.write('cmd.png("conf_%d.png", width=1080, height=1080, dpi=300, ray=1)\n' %i)

cmd = 'pymol -cq friends.py'
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err

cmd = "ffmpeg -framerate 2 -i conf_'%d.png' -plays 0 %s.apng" %args.out
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err
