#!/usr/bin/env python
import argparse, os
import prody
import numpy as np

def protein_align(mobile, static, hetname):
    """
    Align proteins with prody; use a specific chain if necessary and return
    which chain you used, too
    """
    mobile = prody.parsePDB(mobile)
    #see if the desired hetatm residue exists for multiple chains
    chains = np.unique(mobile.select('resname %s' %hetname).getChids())
    if not chains.size:
        print "No chain associated with hetname in ref; attempting to proceed \
        but there may be issues."
    mobile = mobile.select('chain %s ' %chains[0])
    try:
        rmatch, tmatch, seqid, overlap = prody.matchChains(static, mobile,
                pwalign=True, overlap=70)[0]
        moved, transformation = prody.superpose(tmatch, rmatch)
        return mobile,chains[0]
    except:
        print "Failed to match %s, continuing...\n" %pdb
        return prody.parsePDB(pdb),chains[0]

def is_close(coord1, coord2, threshold):
    return ((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 +
            (coord1[2] - coord2[2])**2)**0.5 < threshold

parser = argparse.ArgumentParser(description="Strip water from traj frame based on reference.")

parser.add_argument('-t', '--traj', help='Trajectory frame to strip water from.')
parser.add_argument('-r', '--ref', help='Reference structure consisting of \
        receptor and reference small molecule.')
parser.add_argument('-d', '--distance', type=float, help='Distance threshold for stripping water.')
parser.add_argument('-het', '--hetname', help='Name of small molecule residue to \
        define the stripping distance.')

args = parser.parse_args()

outname = os.path.splitext(args.traj)[0] + '_stripwat.pdb'
#align complex to trajectory frame
static = prody.parsePDB(args.traj)
aligned_complex,chain = protein_align(args.ref, static, args.hetname)
# prody.writePDB('aligned_mobile.pdb', aligned_complex)

#use hetname coords to define threshold for water retention
lig = aligned_complex.select('resname %s and chain %s' %(args.hetname, chain))
# prody.writePDB('aligned_lig.pdb', lig)

kept_res = [] #list of resnums to keep
for traj_atom in static.iterAtoms():
    resnum = traj_atom.getResnum()
    if resnum in kept_res:
        continue
    resname = traj_atom.getResname()
    if resname == 'Cl-' or resname == 'Na+': #remove counterions
        continue 
    if resname != 'WAT' and resname != 'HOH':
        kept_res.append(resnum) #don't touch anything else except waters
    else:
        traj_coord = traj_atom.getCoords()
        for lig_coord in lig.getCoords():
            if is_close(traj_coord, lig_coord, args.distance):
                kept_res.append(traj_atom.getResnum())
                continue
print kept_res
refstring = ''.join(['resnum {} or ' if idx<len(kept_res)-1 else 'resnum {}' for
    idx,item in enumerate(kept_res)])
stripped = static.select(refstring.format(*kept_res))
prody.writePDB(outname, stripped)
