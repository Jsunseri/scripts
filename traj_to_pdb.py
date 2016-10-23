#!/usr/bin/python

import MDAnalysis
import sys,os

u = MDAnalysis.Universe(sys.argv[1],sys.argv[2])
protein = u.select_atoms("protein")
count = 0
base = os.path.basename(sys.argv[1])
for ts in u.trajectory:
    w = MDAnalysis.Writer(base.split('.')[0] + "_traj_" + str(count) + ".pdb", multiframe=True)
    w.write(protein)
    w.close()
    count += 1
