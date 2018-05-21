#!/usr/bin/env python
from rdkit.Chem import AllChem,rdMolAlign
import sys,os

def getmol(file):
    ext = os.path.splitext(file)[-1]
    if ext == ".sdf":
        ref = AllChem.SDMolSupplier(file)
    elif ext == ".pdb":
        ref = AllChem.MolFromPDBFile(file)
    else:
        print "Unsupported ext %s" %ext
        sys.exit()
    return ref

if len(sys.argv[1:]) < 2:
    print "Too few args for alignment.\nUSAGE: shape_align.py [REF] [MOBILE]"

reffile = sys.argv[1]
ref = getmol(reffile)
mobile = getmol(sys.argv[2])
name = os.path.splitext(sys.argv[2])[0]

writer = AllChem.SDWriter('%s_aligned_to_xtal.sdf' %name)
refpymp = AllChem.MMFFGetMoleculeProperties(ref) 
for i,mol in enumerate(mobile):
    prbpymp = AllChem.MMFFGetMoleculeProperties(mol) 
    alignment = rdMolAlign.GetO3A(mol, ref, prbpymp, refpymp)
    score = alignment.Align()
    print '%s,%s,%f' %(name,reffile,score)
    writer.write(mol)

