#!/usr/bin/python

#Generate library files for non-protein residues using antechamber

import os, sys, subprocess, struct
from itertools import izip_longest

#Finds name and file extension of ligand, name of receptor, generates derivative variable names for later


def accumulate(iterable):
	'''
	Makes an accumulator that is used in the fixed-width file parser.
	'''
	total = next(iterable)
	yield total 
	for value in iterable:
		total += value
		yield total 

def make_parser(fieldwidths):
	'''
	Parses a fixed-width formatted file.
	'''
	cuts = tuple(cut for cut in accumulate(abs(fw) for fw in fieldwidths))
	pads = tuple(fw < 0 for fw in fieldwidths)
	flds = tuple(izip_longest(pads, (0,)+cuts,cuts))[:-1]
	parse = lambda line: [line[i:j] for pad, i , j in flds if not pad]
	return parse

#These are the field widths specified for the PDB file format
fieldwidths = (6, 5, 5, 1, 3, 2, 4, 1, 11, 8, 8, 6, 6, 13, 2)

ligand = sys.argv[1]

subprocess.call(["cp", ligand, "LIG_temp.pdb"])
ligand_temp = open(ligand)
os.system("rm -f LIG.pdb")

#Parse the data from the columns
parse = make_parser(fieldwidths)
lig_data = [parse(line) for line in ligand_temp]
ligand_temp.close()
final_ligand = open("LIG.pdb",'w')
count = 0			

for x in range(len(lig_data)):
	if lig_data[x][0]=='ATOM  ' or lig_data[x][0]=='HETATM':
		lig_data[x][4]='LIG' 
		lig_data[x][13] = lig_data[x][13].title()
		lig_data[x][2] = lig_data[x][2].strip()
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
    if lig_data[x][0]=='ATOM  ' or lig_data[x][0]=='HETATM':		      
	lig_data[x][0] = 'HETATM'
        lig_data[x][2] = '{:>{}s}'.format(str(lig_data[x][2]),fieldwidths[2])
	lig_data[x][2] = lig_data[x][2].title()
        final_ligand.write(''.join(lig_data[x]))

final_ligand.write('END\n')
final_ligand.close()	

subprocess.call(["babel", "LIG.pdb", "LIG.sdf"])
#os.system("sed -i 's/%s/LIG/g' LIG.sdf" % (fname))

subprocess.call(["antechamber", "-i", "LIG.sdf", "-fi", "sdf", "-o", "LIG.mol2", "-fo", "mol2", "-c", "gas", "-s", "2"])
os.system("sed -i 's/\<MOL\>/%s/g' %s" % ("LIG", "LIG.mol2"))
subprocess.call(["parmchk", "-i", "LIG.mol2", "-f", "mol2", "-o", "LIG" + ".frcmod"])

leap = open("LIG" + ".in","w")
leap.write('source leaprc.ff14SB\nsource leaprc.gaff\n'+"LIG"+'=loadmol2 '+"LIG.mol2"+'\nloadamberparams '+"LIG"+'.frcmod\n'+'saveoff '+"LIG"+' '+"LIG"+'.lib\nsaveamberparm '+"LIG"+' '+"LIG"+'.prmtop '+"LIG"+'.inpcrd\n'+'quit')
leap.close()
subprocess.call(["tleap","-f","LIG" + ".in"])

os.system("sed -i 's/CL/Cl/g' LIG.pdb")
os.system("rm -f LIG_temp.pdb")
#Prep for MD
#os.system("prepareamber.pl %s %s" % ("LIG.pdb", "LIG"))
