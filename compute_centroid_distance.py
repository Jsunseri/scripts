#!/usr/bin/env python

import pybel
import argparse,os

def get_ftype(fname):
    """ 
    given a molfile, return the ftype for openbabel
    """
    root,ext = os.path.splitext(fname)
    if ext == ".gz":
        ext = os.path.splitext(root)[-1]
    ext = ext.split('.')[-1]
    return ext

def get_centroid(mol):
    """
    given an openbabel mol, compute the centroid
    """
    mol.removeh()
    centroid = [0]*3
    natoms = len(mol.atoms)
    for a in mol.atoms:
        centroid[0] += a.coords[0]
        centroid[1] += a.coords[1]
        centroid[2] += a.coords[2]
    centroid[0] /= natoms
    centroid[1] /= natoms
    centroid[2] /= natoms
    return centroid

def get_len(some_vec):
    """
    given some vector, return its length
    """
    return sum([x*x for x in some_vec])**0.5

parser = argparse.ArgumentParser(description="Compute distance between centroids of ref and other")
parser.add_argument("-r", "--ref", help="Reference molecule whose centroid will be used as the origin for the distance calculation; if more than one conformer, only the first will be used")
parser.add_argument("-m", "--mols", help="Other molecules, whose centroids will be compared with the reference")

args = parser.parse_args()

# read in ref
ext = get_ftype(args.ref)
ref = next(pybel.readfile(ext, args.ref))
origin = get_centroid(ref)

# read in mols
ext = get_ftype(args.mols)
mols = pybel.readfile(ext, args.mols)

# print distances
for mol in mols:
    centroid = get_centroid(mol)
    terms = zip(origin, centroid)
    diffs = [x[1]-x[0] for x in terms]
    print(get_len(centroid))
