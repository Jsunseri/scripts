#!/usr/bin/env python
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
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
    print "Too few args for alignment.\nUSAGE: align_mcs.py [REF] [MOBILE]"

reffile = sys.argv[1]
ref = getmol(reffile)
mobile = getmol(sys.argv[2])
name = os.path.splitext(sys.argv[2])[0]

#for each mol, get MCS with ref and align with obfit.
#unfortunately since obfit takes only a single SMARTS pattern as input, 
#I think I have to use an intermediate temporary flie for each mol in MOBILE
writer = AllChem.SDWriter('%s_aligned_to_xtal.sdf' %name)
for i,mol in enumerate(mobile):
    res = rdFMCS.FindMCS([mol,ref], completeRingsOnly=True)
    core = AllChem.MolFromSmarts(res.smartsString)
    match = mol.GetSubstructMatch(core)   
    refMatch = ref.GetSubstructMatch(core) 
    AllChem.AlignMol(mol, ref, atomMap=zip(match,refMatch)) 
    writer.write(mol)

os.remove(tmpin)
os.remove(tmpout)
