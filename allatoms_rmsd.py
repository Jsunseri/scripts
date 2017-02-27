from pymol import cmd
import sys,os

#load in the structures for RMSD calculation
ref = sys.argv[1]
comp = sys.argv[2]
cmd.load(ref)
cmd.load(comp)

#align them, print output to stdout
ref_name = os.path.splitext(os.path.basename(ref))[0]
comp_name = os.path.splitext(os.path.basename(comp))[0]

rmsd = cmd.align(comp_name, ref_name, cycles=0)[0]
print "%s %s %s\n" % (ref_name, comp_name, rmsd)
