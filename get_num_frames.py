#!/usr/bin/python

import MDAnalysis
import sys,os

assert len(sys.argv[1:]) == 2, "Input: get_num_frames.py <topology> <traj>\n"
assert os.path.splitext(sys.argv[1])[-1] == '.prmtop', "Arg1 should be *.prmtop\n"
assert os.path.splitext(sys.argv[2])[-1] == '.nc', "Arg2 should be *.nc\n"
topology = sys.argv[1]
trajectory = sys.argv[2]

try:
    frames = MDAnalysis.Universe(topology,trajectory)
except:
    print('0')
    sys.exit(1)
print((frames.trajectory.totaltime + frames.trajectory.dt) * 1e-3)
