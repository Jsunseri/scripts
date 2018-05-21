#!/usr/bin/env python
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold
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
    print "Too few args for alignment.\nUSAGE: align_by_scaffold.py [REF] [MOBILE]"

reffile = sys.argv[1]
ref = getmol(reffile)
mobile = getmol(sys.argv[2])
name = os.path.splitext(sys.argv[2])[0]

#for each mol, get MCS of scaffolds and align mols to minimize scaffold RMSD
writer = AllChem.SDWriter('%s_aligned_to_xtal.sdf' %name)
for i,mol in enumerate(mobile):
    refcore = MurckoScaffold.GetScaffoldForMol(ref)
    mobcore = MurckoScaffold.GetScaffoldForMol(mol) 
    res = rdFMCS.FindMCS([refcore,mobcore], completeRingsOnly=True)
    core = AllChem.MolFromSmarts(res.smartsString)
    match = mol.GetSubstructMatch(core)   
    refMatch = ref.GetSubstructMatch(core) 
    AllChem.AlignMol(mol, ref, atomMap=zip(match,refMatch)) 
    writer.write(mol)
