#!/usr/bin/env python

from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
import argparse,os

def parse_file(fname):
    """
    given a filename, use the extension to parse with rdkit; return mol or
    supplier with an indicator of whether it's a supplier or not
    """
    gz = False 
    root,ext = os.path.splitext(fname)
    if ext == ".gz":
        ext = os.path.splitext(root)[-1]
        gz = True
    if gz:
        if ext == ".sdf":
            return (AllChem.SDMolSupplier(fname, strictParsing=False), True)
        else:
            sys.exit("Unsupported extension %s.gz" %ext)
    else:
        if ext == ".sdf" or ext == ".mol":
            return (AllChem.MolFromMolFile(fname, sanitize=False), False)
        elif ext == ".mol2":
            return (AllChem.MolFromMol2File(fname, sanitize=False), False)
        elif ext == ".pdb":
            return (AllChem.MolFromPDBFile(fname, sanitize=False), False)
        else:
            sys.exit("Unsupported extension %s.gz" %ext)

def print_dist(mol, ref):
    """ 
    given two mols, print the distance between their centroids
    """
    origin = ComputeCentroid(ref.GetConformer(0))
    centroid = ComputeCentroid(mol.GetConformer(0))
    print((centroid - origin).length())

parser = argparse.ArgumentParser(description="Compute distance between centroids of ref and other")
parser.add_argument("-r", "--ref", help="Reference molecule whose centroid will be used as the origin for the distance calculation; if more than one conformer, only the first will be used")
parser.add_argument("-m", "--mols", help="Other molecules, whose centroids will be compared with the reference")

args = parser.parse_args()

ref,is_supplier = parse_file(args.ref)
if is_supplier:
    ref = ref.NextMol()

mol,is_supplier = parse_file(args.mols)
if is_supplier:
    for m in mol:
        print_dist(m)
else:
    print_dist(mol)
