#!/usr/bin/python
import os, sys, gzip
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD

ligand = sys.argv[1]
torsions = []
ext = os.path.splitext(ligand)
if ext[-1] == '.sdf':
    f = open(ligand)
elif ext[-1] == '.sdf.gz':
    f = gzip.open(ligand)
else:
    print 'Ligands not provided in sdf format\n'
    sys.exit()

suppl = Chem.ForwardSDMolSupplier(f)
for mol in suppl:
    if mol is None: continue
    torsions.append(rdMD.CalcNumRotatableBonds(mol))
f.close()

plt = sns.distplot(torsions, kde=False, rug=True);
fig = plt.get_figure()
fig.savefig('torsions_dist.pdf', bbox_inches='tight')
print "Mean number of torsions is %f\n" %(sum(torsions)/float(len(torsions)))
