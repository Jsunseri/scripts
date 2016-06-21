#!/usr/bin/python

#Prepare docked complex and ligand sdf for MD. Takes ligand and receptor files as command line arguments, in that order.

#Version 1, July 17, 2015.

import os, sys, subprocess, struct
from itertools import izip_longest

#Finds name and file extension of ligand, name of receptor, generates derivative variable names for later

for arg in sys.argv:
	fname = sys.argv[1].split('.',1)[0]
	base = fname[0:3]
	extension = sys.argv[1].split('.',1)[1]
	receptor = sys.argv[2]
	fixed_receptor = "amber_" + receptor
	mol2 = base + ".mol2"
	sdf = base + ".sdf"
	pdb = base + ".pdb"

#Add in missing hydrogens with babel, generate mol2 file with 
#antechamber, substitute ligand name into file, generate library files
subprocess.call(["babel", "-h", fname + '.' + extension, sdf])
subprocess.call(["babel", "-h", fname + '.' + extension, pdb])
subprocess.call(["antechamber", "-i", sdf, "-fi", "sdf", "-o", mol2, "-fo", "mol2", "-c", "bcc", "-s", "2"])
os.system("sed -i 's/\<MOL\>/%s/g' %s" % (base, mol2))
subprocess.call(["parmchk", "-i", mol2, "-f", "mol2", "-o", base + ".frcmod"])
os.system("pdb4amber -y -p -i %s -o %s" % (receptor, fixed_receptor))

leap = open(base + ".in","w")
leap.write('source leaprc.ff14SB.redq\nsource leaprc.gaff\n'+base+'=loadmol2 '+mol2+'\nloadamberparams '+base+'.frcmod\n'+'saveoff '+base+' '+base+'.lib\nsaveamberparm '+base+' '+base+'.prmtop '+base+'.inpcrd\n'+'quit')
leap.close()
subprocess.call(["tleap","-f",base + ".in"])

#Generate pdb of receptor-ligand complex
f = open(fixed_receptor,'r')
g = open(base + receptor,'w')

for line in f:
	if not line.startswith('END'):
		g.write(line)
g.write('TER\n')

f.close()
g.close()

os.system("sed -e 's/LIG/%s/g' -e 's/ATOM[ ][ ]/HETATM/g' -e '/HETATM/!d' <%s >temp_lig.pdb" % (base, pdb))
os.system("sed '/ATOM/!d' <%s >temp_rec.pdb" % ("amber_" + receptor))

#Set up functions to analyze and generate fixed-width files

def accumulate(iterable):
	total = next(iterable)
	yield total 
	for value in iterable:
		total += value
		yield total 

def make_parser(fieldwidths):
	cuts = tuple(cut for cut in accumulate(abs(fw) for fw in fieldwidths))
	pads = tuple(fw < 0 for fw in fieldwidths)
	flds = tuple(izip_longest(pads, (0,)+cuts,cuts))[:-1]
	parse = lambda line: [line[i:j] for pad, i , j in flds if not pad]
	return parse

fieldwidths = (6, 5, 5, 1, 3, 2, 4, 1, 11, 8, 8, 6, 6, 13, 2)

f = open("temp_rec.pdb","rw")
g = open("temp_lig.pdb","rw")

parse = make_parser(fieldwidths)
rec_data = [parse(line) for line in f]
lig_data = [parse(line) for line in g]

#Update the atom numbers, residue numbers, and atom names for the ligand in the complex PDB file

for x in range(len(lig_data)):
	lig_data[x][1] = int(rec_data[-1][1]) + int(lig_data[x][1])
	lig_data[x][6] = int(rec_data[-1][6]) + 1
	lig_data[x][2] = lig_data[x][2].strip()
	lig_data[x][2] = lig_data[x][2].title()
	lig_data[x][2] = ''.join([i for i in lig_data[x][2] if not i.isdigit()])
	
occurrences = {}
for index, atom in enumerate(lig_data[x][2] for x in range(len(lig_data))):
	if atom not in occurrences:
		occurrences[atom] = [index,1]
		lig_data[index][2] += str(occurrences[atom][1])
	else:
		occurrences[atom][1] += 1
		lig_data[index][2] += str(occurrences[atom][1])

for x in range(len(lig_data)):
	lig_data[x][1] = '{:>{}s}'.format(str(lig_data[x][1]),fieldwidths[1])
	lig_data[x][6] = '{:>{}s}'.format(str(lig_data[x][6]),fieldwidths[6])
	lig_data[x][2] = '{:>{}s}'.format(str(lig_data[x][2]),fieldwidths[2])
			
with open(base + receptor,'a') as doc:
	for x in range(len(lig_data)):
		doc.write(''.join(lig_data[x]))
	doc.write('END\n')

os.system("rm -f temp*")

#Prep for MD
os.system("prepareamber.pl %s %s" % (base + receptor, base))
