#!/usr/bin/env python
import argparse, os
import pdb_util as util
#openbabel segfaults on these examples with water boxes; rdkit seems impossibly
#slow; and biopython spuriously complains about repeated residues for the waters

def get_moldata(line):
    mol_data = {}
    for i,field in enumerate(util.pdb_fieldnames):
        if i in util.pdb_floatfields:
            field_info = float(line[i]) if line[i].strip() else line[i]
        elif i in util.pdb_intfields:
            field_info = int(line[i]) if line[i].strip() else line[i]
        else:
            field_info = line[i] 
        mol_data[field] = field_info
    return mol_data

def write_line(rec_contents, lig_coords, parse, distance):
    if "TER" in line:
        return False
    closest_dist = 1000
    for coord in lig_coords:
        dist = ((coord['x'] - rec_contents['x'])**2 + (coord['y'] -
            rec_contents['y'])**2 + (coord['z'] - rec_contents['z'])**2)**0.5
        if dist < closest_dist: closest_dist = dist
    return closest_dist <= distance

parser = argparse.ArgumentParser(description="Strip water from traj frame based on reference.")

parser.add_argument('-t', '--traj', help='Trajectory frame to strip water from.')
parser.add_argument('-r', '--ref', help='Reference to define the stripping distance.')
parser.add_argument('-d', '--distance', type=float, help='Distance threshold for stripping water.')

args = parser.parse_args()
parse = util.make_parser(util.pdb_fieldwidths)

rec = open(args.traj, 'r')
lig = open(args.ref, 'r')
outname = os.path.splitext(args.traj)[0] + '_stripwat.pdb'
out = open(outname, 'w')
lig_coords = []
for line in lig:
    if not "ATOM" in line and not "HETATM" in line:
        continue
    lig_contents = get_moldata(parse(line))
    lig_coords.append({'x': lig_contents['x'], 'y': lig_contents['y'], 'z':
        lig_contents['z']})

kept_waters = [] #list of water resnums we've decided to keep
for line in rec: 
    if "Cl-" in line or "Na+" in line: #remove counterions
        continue
    if not "WAT" in line: #don't touch anything but waters
        out.write(line)
    else:
        contents = get_moldata(parse(line))
        if contents['resnum'] in kept_waters:
            out.write(line)
        elif write_line(contents, lig_coords, parse, args.distance):
            kept_waters.append(contents['resnum'])
            out.write(line)

rec.close()
lig.close()
out.close()
