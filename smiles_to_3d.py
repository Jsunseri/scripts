#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem
import sys,os

# Takes an input file of Smiles strings and outputs an sdf with a reasonable 3d
# conformation for each molecule suitable for something like flexible-ligand
# docking. Takes the file name as the first argument, and optionally a second
# argument that specifies an output file name (without the extension).
# Otherwise it infers the filename from the input filename.

fname = sys.argv[1]
if len(sys.argv) > 2:
    outname = sys.argv[2]
else:
    outpath,outbase = os.path.split(fname)
    outname = os.path.join(outpath, outbase.split('.')[0])

w = Chem.SDWriter(outname+'.sdf')
count = 0

with open(fname,'r') as f:
    for line in f:
        m = Chem.MolFromSmiles(line)
        if m is None:
            print "NoneMol returned for %s, continuing\n" %(line)
            continue
        m2 = Chem.AddHs(m)
        result = AllChem.EmbedMolecule(m2)
        if result != 0:
            result = AllChem.EmbedMolecule(m2,useRandomCoords=True)
            if result != 0:
                print "Failed to embed $s, continuing\n" %(line)
                continue
        AllChem.UFFOptimizeMolecule(m2)
        m2.SetProp("_Name", "mol"+str(count))
        w.write(m2)
        count += 1

w.close()
