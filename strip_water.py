#!/usr/bin/env python
import argparse, os

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

parser = argparse.ArgumentParser(description="Strip water from traj frame based on reference.")

parser.add_argument('-t', '--traj', help='Trajectory frame to strip water from.')
parser.add_argument('-r', '--ref', help='Reference structure consisting of'
        ' receptor and reference small molecule.')
parser.add_argument('-l', '--lig', help='Reference ligand if it\'s not in ref.')
parser.add_argument('-d', '--distance', type=float, help='Distance threshold for stripping water.')
parser.add_argument('-het', '--hetname', help='Name of small molecule residue to'
        ' define the stripping distance; not necessary if you provided a lig.')

args = parser.parse_args()

outname = os.path.splitext(args.traj)[0] + '_stripwat.pdb'
#align complex to trajectory frame
mobile = os.path.splitext(os.path.basename(args.traj))[0]
static = os.path.splitext(os.path.basename(args.ref))[0]
pymol.cmd.load(args.traj, mobile)
pymol.cmd.load(args.ref, static)
if args.lig:
    pymol.cmd.load(args.lig, 'lig')
else:
    assert args.hetname, "If you don't pass a lig you gotta pass a hetname"
    pymol.cmd.select('lig', '%s and resn %s' %(static, args.hetname))
pymol.cmd.select('nowat', '%s and polymer' %mobile)
pymol.cmd.align('nowat', static, object='align')

pymol.cmd.select('strip', '(%s and (not resn WAT and not resn Cl- and not resn Na+)) or (%s and resn WAT within %f of lig)' %(mobile, mobile, args.distance))
pymol.cmd.save(outname, 'strip')
