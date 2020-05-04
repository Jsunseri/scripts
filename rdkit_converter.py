#!/usr/bin/python3
import sys,string,argparse
from rdkit.Chem import AllChem as Chem
import sys, os, gzip

args = sys.argv[1:]
if(len(args) < 2):
    print("Need input and output")
    sys.exit(-1)
    
input = args[0]
output = args[1]
itype = ''
otype = ''

# deal with input; right now just deal with smi or mol2 input
split = os.path.splitext(input)
if split[1] == ".smi":
    itype = "smi"
    ifile = open(input)
elif split[1] == ".mol2":
    itype = "mol2"
    mol = Chem.MolFromMol2File(input)
elif split[1] == ".pdb":
    itype = "pdb"
    mol = Chem.MolFromPDBFile(input)
else:
    assert 0, "Only .smi or .mol2 or .pdb input currently supported due to laziness"

split = os.path.splitext(output)
if split[1] == '.gz':
    outf=gzip.open(output,'w+')
    output = split[0] #strip .gz
else:
    outf = open(output,'w+')
 
if os.path.splitext(output)[1] == '.pdb':
    owriter = Chem.PDBWriter(outf)
elif os.path.splitext(output)[1] == '.sdf':
    owriter = Chem.SDWriter(outf)
elif os.path.splitext(output)[1] == ".smi":
    owriter = Chem.SmilesWriter(outf)
else:
    assert 0, "Only .pdb or .sdf output currently supported"
    
if owriter is None:
    print("Could not open ".output)
    sys.exit(-1)

if itype == "smi":
    for line in ifile:
        toks = line.split()
        smi = toks[0]
        name = ' '.join(toks[1:])    
        
        pieces = smi.split('.')
        if len(pieces) > 1:
            smi = max(pieces, key=len) #take largest component by length
            print("Taking largest component: %s\t%s" % (smi,name))
            
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
                mol = Chem.AddHs(mol)
                mol.SetProp("_Name",name);
                conf = Chem.EmbedMultipleConfs(mol, 1,randomSeed=0)
                Chem.UFFOptimizeMolecule(mol,confId=0)
                mol = Chem.RemoveHs(mol)
                owriter.write(mol,conf)
            except (KeyboardInterrupt, SystemExit):
                raise                
            except Exception as e:
                print("Exception",e)
        else:
            print("ERROR:",smi)
elif itype == "mol2" or itype == "pdb":
    owriter.write(mol)

owriter.close()
outf.close()
