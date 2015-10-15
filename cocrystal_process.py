#!/usr/bin/python

#Decompose PDB file containing information for a protein and multiple ligands/ions/solvent molecules
#into multiple properly formatted files to begin an MD run. Takes the full PDB file containing
#ligands as a command line argument. 

#Version 1, July 22, 2015.

import sys, os, mechanize
from itertools import izip_longest

def retrieve_bindinfo(fname):
	'''
	Get the text of the PDB summary page for the structure specified by 
	fname, which should be the PDBID for the protein of interest.
	'''
	br = mechanize.Browser()
	br.set_handle_robots(False)
	br.set_handle_refresh(False)
	response = br.open('http://www.rcsb.org/pdb/explore/explore.do?structureId=' + fname)
	return response.read()

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

pdb = sys.argv[1]
fname = pdb.split('.',1)[0]
	
#Separate the ATOM and HETATM records
with open(pdb, 'r') as infile, open("rec_temp.pdb",'a+') as receptor_temp, open("lig_temp.pdb",'a+') as ligand_temp:
	for line in infile:
		if line.startswith("ATOM"):
			receptor_temp.write(line)
		if line.startswith("HETATM") and "HOH" not in line:
			ligand_temp.write(line)

ligand_temp = open("lig_temp.pdb",'a+')
receptor_temp = open("rec_temp.pdb",'a+')

#Parse the data from the columns
parse = make_parser(fieldwidths)
rec_data = [parse(line) for line in receptor_temp]
lig_data = [parse(line) for line in ligand_temp]
receptor_temp.close()
ligand_temp.close()

#Identify the unique ligands in the file
fragments = {}
for x in range(len(lig_data)):
	if lig_data[x][4] not in fragments:
		fragments[lig_data[x][4]] = lig_data[x][5]

#If multiple ligands are found, check the PDB website to see which have
#binding affinity data available
if len(fragments) > 1: 		
	choice = None
	website_full = retrieve_bindinfo(fname).splitlines()
	website_subset = []

	copy = False
	for line in website_full:
		if line.strip() == "Binding Affinity (Sequence Identity %)":
			copy = True
		elif line.strip() == "<b>External Domain Annotations</b>":
			copy = False
		elif copy:
			website_subset.append(line)

	crystal_ligs = []			
	for element in website_subset:
		if "Download ligand SDF file for " in element:
			pared_line = element.strip()
			crystal_ligs.append(pared_line.split()[5][0:3])

#Prompts the user to select the desired ligand if the PDB doesn't provide
#binding affinity data or if it has data for multiple ligands in the file			
	while choice is None:
		if len(crystal_ligs) == 0 or len(crystal_ligs) > 1:
			print "Multiple ligands found; please choose one" 
			print [key for key in fragments.keys()]
			choice = raw_input()
			if choice not in fragments:
				print "Invalid selection"
				choice = None
		elif len(crystal_ligs) == 1:
			choice = crystal_ligs[0]
else:
	choice = fragments.keys()[0]

#Creates a file LIG.pdb containing the structure of the chosen ligand.
final_ligand = open("LIG.pdb",'w')
count = 0			
for x in range(len(lig_data)):
	if lig_data[x][4]==choice and lig_data[x][5]==fragments[choice]:
		lig_data[x][4]='LIG'
		lig_data[x][13].title()
		lig_data[x][2] = lig_data[x][13].strip() + str(count)
		lig_data[x][2] = '{:>{}s}'.format(str(lig_data[x][2]),fieldwidths[2])
		final_ligand.write(''.join(lig_data[x]))
		count += 1
final_ligand.write('END\n')
final_ligand.close()	

#Creates a file receptor.pdb containing the receptor without the ligands.
final_rec = open("receptor.pdb",'w')		
for y in range(len(rec_data)):
	if rec_data[y][5]==fragments[choice]:
		final_rec.write(''.join(rec_data[y]))

final_rec.close()
os.system("rm -f *_temp.*")
