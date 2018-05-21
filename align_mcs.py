#!/usr/bin/env python
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import sys,os
import filecmp
try:
    from plumbum.cmd import obfit
except ImportError:
    raise ImportError('Check that obfit is on your path')

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
    print "Too few args for alignment.\nUSAGE: alignmols.py [REF] [MOBILE]"

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
    smarts = res.smartsString
    tmpin = 'intmp.sdf'
    tmpwriter = AllChem.SDWriter(tmpin)
    tmpwriter.write(mol)
    out = obfit[smarts, reffile, tmpin]()
    if len(out) != 0:
        tmpout = 'outtmp.sdf'
        with open(tmpout, 'w') as outf:
            for line in out:
                outf.write(line)
        if filecmp.cmp(tmpin, tmpout, shallow=False):
            print 'Failed to align mol %d\n' %i
            continue
        else:
            tmpmol = AllChem.SDMolSupplier(tmpout)[0]
            writer.write(tmpmol)

os.remove(tmpin)
os.remove(tmpout)
