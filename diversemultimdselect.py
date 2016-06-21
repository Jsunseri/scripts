#!/usr/bin/python

import MDAnalysis, argparse
from MDAnalysis.analysis.align import rmsd
import MDAnalysis.analysis
import os
from os.path import splitext, basename

#This finds the most distinct set of frames.  Does not perform alignment, just rmsd.

def compute_cluster_sizes(current, model, selection):
	"""Given a set of frames (current), assign each from in
	the trajectory to its closest current frame to get cluster sizes"""
	cnts = [0]*len(current)
	modelnumber = 0
	for eachmodel in model:
		for ts in eachmodel.trajectory:
			frame = ts.frame
			minr = float('infinity')
			closest = 0
			for (pos,(everymodel, refframe, refcoords)) in enumerate(current):
				r = rmsd(refcoords, selection[modelnumber].coordinates())			
				if r < minr:
					minr = r
					closest = pos
			cnts[closest] += 1
		modelnumber = modelnumber + 1
	return cnts
		
def add_next_farthest(current, model, selection):
	"""Given a list of frames (frame number, coordinates) current
		identify the frame that has the maximum minimum distance
		from current and add it to current"""
	maxr = 0
	maxframe = 0
	maxcoords = None
	modelnumber = 0
	maxmodel = 0
	for eachmodel in model:
		for ts in eachmodel.trajectory:
			frame = ts.frame
			minr = float('infinity')
			for (everymodel, refframe, refcoords) in current:
				r = rmsd(refcoords, selection[modelnumber].coordinates())			
				if r < minr:
					minr = r
			if minr > maxr:
				maxframe = frame
				maxr = minr
				maxcoords = selection[modelnumber].coordinates()
				maxmodel = modelnumber
		modelnumber = modelnumber + 1
	
	current.append( (maxmodel, maxframe, maxcoords) )

parser = argparse.ArgumentParser(description='Identify most distinct frames of md trajectory\nIMPORTANT: assumes an aligned input')
parser.add_argument("--topologies", nargs='+')
parser.add_argument("--trajectories", nargs='+')
parser.add_argument("--size",type=int)
parser.add_argument("--selections", nargs='+', default="backbone",required=False)
args = parser.parse_args()

model = []
selection = []
modelnumber = 0
for topology, trajectory, selections in zip(args.topologies,args.trajectories,args.selections):
	model.append(MDAnalysis.Universe(topology,trajectory))
	selection.append(model[modelnumber].selectAtoms(args.selections[modelnumber]))
	modelnumber = modelnumber + 1

# current is a list of tuples tracking model number, frame number, and coordinates of the selected atoms for the cluster centers
current = [(0, 1, selection[0].coordinates())]

for i in xrange(args.size-1):
	add_next_farthest(current, model, selection)

name = splitext(basename(args.topologies[0]))[0][0:3]
for (eachmodel,frame,coords) in current:
	model[eachmodel].trajectory[frame-1] #has side effect of setting current frame
	fname = "%s_model%d_%d.pdb" % (name, eachmodel, frame)
	model[eachmodel].selectAtoms('all').write(fname)
	os.system("sed -i '/REMARK/d' %s" % fname) #remove remarks with binary characters

#print out frame numbers and their respective cnts
cnts = compute_cluster_sizes(current,model,selection)
for ((eachmodel,frame,coords),cnt) in zip(current,cnts):
	print 'model %d, %d:%d' % (eachmodel,frame,cnt),
print ''
