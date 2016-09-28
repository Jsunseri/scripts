#!/usr/bin/python

# Creates animated GIF of systematic sampling of subset of ligand parameters.
# version 1
# Needs some work, d00d.

import sys, os, math, random
from rdkit.Chem import AllChem
import subprocess as sp

def make_directory(base):
	""" 
	Make directory with error checking that 
	avoids overwriting an existing directory.
	"""
		
	i = 0
	while 1:
	    try:
	        if i == 0:
	            dirname = base
	        else:
	            dirname = base + '_' + str(i)
	        os.mkdir(dirname)
	        break
	    except OSError:
	        if not os.path.isdir(dirname):
	            raise
	        i += 1
	        pass 
	return dirname

def initialize_conf(old_conf, new_conf):
	"""
	Initializes new conformer with an old one.
	"""
	for atom in range(0, old_conf.GetNumAtoms()):
		old_position = old_conf.GetAtomPosition(atom)
		new_conf.SetAtomPosition(atom, old_position)
	return

# Generate conformers that sample a subset of parameters
suppl = AllChem.SDMolSupplier(sys.argv[1])
mol = suppl[0]
params = {}
param_idx = {}
all_params = ['bond','angle','dihe']
# We'll use these later - the params dictionary will link parameter types to
# lists of the relevant atoms in the molecule of interest and the param_idx
# dictionary will give a lookup for the index that corresponds to that
# parameter type in the collective data structures
for idx,param in enumerate(all_params):
	params[param] = []
	param_idx[param] = idx

param_range = [(1.0,2.0),(40,320),(-180,180)]
param_increment = [0.25,10,10]

num_steps_needed = [(elem[1]-elem[0])/param_increment[i] for i,elem in enumerate(param_range)]
num_steps_needed = [x*y for x,y in zip(num_steps_needed,[2,2,1])]
max_steps = max(num_steps_needed)

# This is the multiplicative factor for the number of steps each parameter type will actually take
num_steps = [math.floor(max_steps/nsteps) for nsteps in num_steps_needed]

adjacency_list = [[] for i in range(mol.GetNumAtoms())]

# Populate adjacency list
for atom in mol.GetAtoms():
	for bond in atom.GetBonds():
		atom_idx = atom.GetIdx()
		if atom_idx == bond.GetBeginAtomIdx():
			adjacency_list[atom_idx].append(bond.GetEndAtomIdx())

# TODO: exclude rings except at end points
for idx,sublist in enumerate(adjacency_list):
	for child in sublist:
		if not mol.GetAtomWithIdx(idx).IsInRing() or not mol.GetAtomWithIdx(child).IsInRing():
			params['bond'].append((idx,child))
			for gchild in adjacency_list[child]:
				if not mol.GetAtomWithIdx(gchild).IsInRing():
					params['angle'].append((idx,child,gchild))
					for ggchild in adjacency_list[gchild]:
						params['dihe'].append((idx,child,gchild,ggchild))

chosen_bond = params['bond'][random.randint(0,len(params['bond'])-1)]
chosen_angle = params['angle'][random.randint(0,len(params['angle'])-1)]
chosen_dihe = params['dihe'][random.randint(0,len(params['dihe'])-1)]

bond_mode = 1
angle_mode = 1

dirname = make_directory('ligand_GIF')
writer = AllChem.SDWriter(dirname+'/conformers.sdf')
writer.write(mol,0)

# TODO: clean this up, geez Louise
for i in range(0,max_steps):
	old_conf = mol.GetConformer(i)
	new_conf = AllChem.Conformer(mol.GetNumAtoms())
	initialize_conf(old_conf, new_conf)
	new_confId = mol.AddConformer(new_conf,True)
        current_conf = mol.GetConformer(new_confId)
	if i < num_steps[param_idx['bond']] * num_steps_needed[param_idx['bond']]:
		old_length = AllChem.GetBondLength(old_conf,chosen_bond[0],chosen_bond[1])
		new_length = old_length + param_increment[param_idx['bond']] * bond_mode
		if new_length > param_range[param_idx['bond']][0] and new_length < param_range[param_idx['bond']][1]:
			AllChem.SetBondLength(current_conf,chosen_bond[0],chosen_bond[1],new_length)
		else:
			bond_mode = bond_mode * -1
			new_length = old_length + param_increment[param_idx['bond']] * bond_mode
			AllChem.SetBondLength(current_conf,chosen_bond[0],chosen_bond[1],new_length)
	if i < num_steps[param_idx['angle']] * num_steps_needed[param_idx['angle']]:
		old_theta = AllChem.GetAngleDeg(new_conf,chosen_angle[0],chosen_angle[1],chosen_angle[2])
		new_theta = old_theta + param_increment[param_idx['angle']] * angle_mode
		if new_theta > param_range[param_idx['angle']][0] and new_theta < param_range[param_idx['angle']][1]:
			AllChem.SetAngleDeg(current_conf,chosen_angle[0],chosen_angle[1],chosen_angle[2],new_theta)
		else:
			angle_mode = angle_mode * -1
			new_theta = old_theta + param_increment[param_idx['angle']] * angle_mode
			AllChem.SetAngleDeg(current_conf,chosen_angle[0],chosen_angle[1],chosen_angle[2],new_theta)
	if i < num_steps[param_idx['dihe']] * num_steps_needed[param_idx['dihe']]:
		old_tor = AllChem.GetDihedralDeg(new_conf,chosen_dihe[0],chosen_dihe[1],chosen_dihe[2],chosen_dihe[3])
		new_tor = old_tor + param_increment[param_idx['dihe']]
		AllChem.SetDihedralDeg(current_conf,chosen_dihe[0],chosen_dihe[1],chosen_dihe[2],chosen_dihe[3],new_tor)
	writer.write(mol,i+1)

os.chdir(dirname)
with open('friends.py','w') as file:
    file.write('from pymol import cmd\n\n')
    file.write('cmd.load("conformers.sdf")\n')
    file.write('cmd.bg_color("white")\n')
    file.write('cmd.show_as("sticks")\n')
    #file.write('cmd.set("ray_trace_frames","1")\n')
    for i in range(1,max_steps+1):
        file.write('cmd.set("state",%d)\n' %i)
        file.write('cmd.ray()\n')
        file.write('cmd.png("conf_%d.png")\n' %i)

cmd = 'pymol -cq friends.py'
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err

cmd = 'convert -delay 20 -loop 0 *.png ligand_friend.gif'
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err
