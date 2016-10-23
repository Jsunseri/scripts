#!/usr/bin/python
import sys,os
from rdkit import Chem

if os.path.isfile(sys.argv[1]):
    base,ext = sys.argv[1].split('.')
    with open(sys.argv[1],'r') as f:
        with open(base + '_canon.' + ext,'w') as g:
            for line in f:
                g.write(Chem.MolToSmiles(Chem.MolFromSmiles(line)) + '\n')
else:
    print Chem.MolToSmiles(Chem.MolFromSmiles(sys.argv[1]))
