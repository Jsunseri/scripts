#!/usr/bin/python
from Bio.PDB.Polypeptide import is_aa
from itertools import izip_longest

pdb_fieldwidths = (6, 5, 5, 1, 3, 2, 4, 1, 11, 8, 8, 6, 6, 13, 2)

def accumulate(iterable):
	total = next(iterable)
	yield total 
	for value in iterable:
		total += value
		yield total 

def make_parser(pdb_fieldwidths):
	cuts = tuple(cut for cut in accumulate(abs(fw) for fw in pdb_fieldwidths))
	pads = tuple(fw < 0 for fw in pdb_fieldwidths)
	flds = tuple(izip_longest(pads, (0,)+cuts,cuts))[:-1]
	parse = lambda line: [line[i:j] for pad, i , j in flds if not pad]
	return parse

def parsepdb(pdb):
    parse = make_parser(pdb_fieldwidths)
    f = open(pdb, 'r')
    mol_data = [parse(line) for line in f if line.startswith('HETATM') or line.startswith('ATOM')]
    f.close()
    return mol_data

def reorder_residues(mol_data):
    resmap = {}
    for x in range(len(mol_data)):
        resnum = mol_data[x][6].strip()
        if resnum not in resmap.keys():
            resmap[resnum] = [x]
        else:
            resmap[resnum].append(x)
    x = 1
    new_indices = []
    for num in range(1,len(resmap.keys())+1):
        for atom in resmap[str(num)]:
            mol_data[atom][1] = '{:>{}s}'.format(str(x),pdb_fieldwidths[1])
            new_indices.append(atom)
            x += 1
    new_mol_data = [mol_data[i] for i in new_indices]
    return new_mol_data

def renumber_atoms(mol_data):
    for x in range(len(mol_data)):
        mol_data[x][1] = x + 1
        mol_data[x][1] = '{:>{}s}'.format(str(mol_data[x][1]),pdb_fieldwidths[1])
    return mol_data

def fix_unlabeled_chains(mol_data):
    res_seen = {}
    greedy_mods = {}
    for x in range(len(mol_data)):
        res_id = mol_data[x][6]
        chain = mol_data[x][5].strip()
        if not chain:
            if res_id in res_seen:
                chain = res_seen[res_id]
                if chain in greedy_modes.keys():
                    greedy_mods[chain].append(x)
            else:
                if greedy_mods.keys():
                    # Technically this *could* be a problem if there were a lot
                    # of atoms with unmarked chains
                    chain = chr(ord(greedy_mods.keys()[-1]) + 1)
                else:
                    # Chosen because this never appeared as a chain identifier
                    # for any of the files on my computer
                    chain = 'f'
                res_seen[res_id] = chain
                greedy_mods[chain] = [x]
            mol_data[x][5] = '{:>{}s}'.format(chain,fieldwidths[5])
        else:
            res_seen[res_id] = chain
            if chain in greedy_mods.keys():
                new_chain = chr(ord(greedy_mods.keys()[-1]) + 1)
                greedy_mods[new_chain] = []
                for atom in greedy_mods[chain]:
                    res_seen[mol_data[atoms][6]] = new_chain
                    greedy_mods[new_chain].append(atom)
    outfile = ''
    # If a fix was required, you need to output a new complex file with the
    # updated chain information; returning the file name (which otherwise will
    # be an empty string) lets the caller know what you did
    if greedy_mods:
        outfile = 'complex_fixedchains.pdb'
        writepdb(mol_data,outfile)

    return mol_data,outfile

def rename_atoms(mol_data):
    for x in range(len(mol_data)):
    	mol_data[x][2] = mol_data[x][2].strip()
    	mol_data[x][2] = mol_data[x][2].title()
    	mol_data[x][2] = ''.join([i for i in mol_data[x][2] if not i.isdigit()])
    	
    occurrences = {}
    for index, atom in enumerate(mol_data[x][2] for x in range(len(mol_data))):
    	if atom not in occurrences:
    		occurrences[atom] = [index,1]
    		mol_data[index][2] += str(occurrences[atom][1])
    	else:
    		occurrences[atom][1] += 1
    		mol_data[index][2] += str(occurrences[atom][1])
        mol_data[index][2] = '{:>{}s}'.format(str(mol_data[index][2]),fieldwidths[2])
    return mol_data

def not_prot(res,record_type):
    if not any([is_aa(res), res=='HIE', res=='HID', res=='CYX', res=='CYM',
        res=='WAT', res=='HOH']) or record_type == 'HETATM' and not any([res=='WAT', res=='HOH']):
        return True
    else:
        return False

def get_nonprot(mol_data):
    nonprot = {}
    for x in range(len(mol_data)):
        res = mol_data[x][4]
        resnum = mol_data[x][6]
        if not_prot(res, mol_data[x][0]): 
            if resnum not in nonprot.keys():
                nonprot[resnum] = [res,1]
            else:
                nonprot[resnum][1] += 1
    resnums = nonprot.keys()
    resnames = [i[0] for i in nonprot.values()]
    natoms = [i[1] for i in nonprot.values()]
    return resnums,resnames,natoms

def get_chains(mol_data):
    rec_chains  = []
    lig_chains = {}
    for x in range(len(mol_data)):
        chain_id = mol_data[x][5]
        if not_prot(mol_data[x][4], mol_data[x][0]):
            resnum = mol_data[x][6]
            if resnum not in lig_chains.keys():
                lig_chains[resnum] = chain_id
        else:
            if chain_id not in rec_chains:
                rec_chains.append(chain_id)
    return rec_chains,lig_chains

def writepdb(mol_data, fname):
    with open(fname, 'w') as f:
        for i in range(len(mol_data)):
            mol_data[x][13] = mol_data[x][13].title()
            mol_data[x][13] = '{:>{}s}'.format(str(mol_data[x][13]),fieldwidths[13])
            f.write(''.join(mol_data[i]))
        f.write('END\n')
    return
