#!/usr/bin/python
from itertools import izip_longest

#indices 2 6 10 16 are padding fields associated with 
#columns not assigned to anything by the spec
pdb_fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2)
pdb_fieldnames = ('recordname', 'atomnum', 'atomname', 'altloc',
        'resname', 'chainid', 'resnum', 'rescode', 'x', 'y', 'z', 'occupancy',
        'beta', 'element', 'charge')

#amino acid residues that leap should recognize with a standard protein force
#field, plus water
standard_res = set([ALA, GLY, SER, THR, LEU, ILE, VAL, ASN, GLN, ARG, HID, HIE,
    HIP, TRP, PHE, TYR, GLU, ASP, LYS, LYN, PRO, CYS, CYX, MET, ASH, GLH, CYM,
    HYP, NALA, NGLY, NSER, NTHR, NLEU, NILE, NVAL, NASN, NGLN, NARG, NHID, 
    NHIE, NHIP, NTRP, NPHE, NTYR, NGLU, NASP, NLYS, NPRO, NCYS, NCYX, NMET,
    NME, NHE, ACE, CALA, CGLY, CSER, CTHR, CLEU, CILE, CVAL, CASN, CGLN, CARG, 
    CHID, CHIE, CHIP, CTRP, CPHE, CTYR, CGLU, CASP, CLYS, CPRO, CCYS, CCYX, 
    CMET, CHYP, WAT, HOH])

#ions and cofactors that have parameters at the AMBER parameter database
#http://research.bmh.manchester.ac.uk/bryce/amber
#and a selection of other common molecules like SO4 that also have widely
#available parameters that aren't distributed with AMBER
common_ions_and_cofactors = set([GDP, GTP, ADP, ATP, FMN, FAD, NAD, NADH, NAH,
    NDP, NPD, NPH, ARP, HEM, HEME, MG, CO6, CAL, MNG, SO4, PO4])

#modified amino acid residues that have parameters at the AMBER parameter
#database http://research.bmh.manchester.ac.uk/bryce/amber
modified_residues = set([ZA, ZC, ZD, ZE, ZF, ZG, ZHD, ZHE, ZHP, ZI, ZK, ZL, ZM, 
    ZN, ZP, ZQ, ZR, ZS, ZT, ZV, ZW, ZY, M3L, CC4, CH4, K3M, TFL, HFL, NOR, ORN,
    HEP, H1D, H2D, H1E, H2E, S1P, S2P, T1P, T2P, Y1P, Y2P, SEP, THP, TYP])

def accumulate(iterable):
    '''
    Generate running sum from an iterator
    '''	
    total = next(iterable)
    yield total 
    for value in iterable:
    	total += value
    	yield total 

def make_parser(pdb_fieldwidths):
    '''
    Return a fixed-width file parser
    '''
    cuts = tuple(cut for cut in accumulate(abs(fw) for fw in pdb_fieldwidths))
    pads = tuple(fw < 0 for fw in pdb_fieldwidths)
    flds = tuple(izip_longest(pads, (0,)+cuts, cuts))[:-1]
    parse = lambda line: [line[i:j].rstrip('\n') for pad, i, j in flds if not pad]
    return parse

def parsepdb(pdb):
    '''
    Parse a PDB file and stores the ATOM and HETATM records in a dictionary
    '''
    parse = make_parser(pdb_fieldwidths)
    f = open(pdb, 'r')
    mol_data_list = [parse(line) for line in f if line.startswith('HETATM') or line.startswith('ATOM')]
    f.close()
    mol_data = {}
    for i,field in enumerate(pdb_fieldnames):
        if (i == 1 or i == 6 or i == 8 or i == 9 or i == 10 or i == 11 or i ==
        12 or i == 14):
            fieldlist = [float(line[i]) for line in mol_data_list]
        else:
            fieldlist = [line[i] for line in mol_data_list]
        mol_data[field] = fieldlist
    if not mol_data['element']:
        mol_data = set_element(mol_data)
    return mol_data

def group_by_residue(mol_data):
    '''
    Rearrange atoms in a file so that atoms in the same residue are
    contiguous
    '''
    resmap = {}
    for old_idx in range(len(mol_data.atomnum)):
        resnum = mol_data.resnum[old_idx]
        if resnum not in resmap.keys():
            resmap[resnum] = [old_idx]
        else:
            resmap[resnum].append(old_idx)
    atnum = 1
    new_indices = [None] * len(mol_data.atomnum)
    for res in resmap.keys():
        for old_idx in resmap[res]:
            mol_data.atomnum[old_idx] = atnum
            new_indices.append(atnum-1)
            atnum += 1
    new_mol_data = {}
    for key in mol_data.keys():
        new_mol_data[key] = [mol_data[key][i] for i in new_indices]
    return new_mol_data

def renumber_atoms(mol_data, start_val=1):
    '''
    Renumber atoms so they start at 1
    '''
    for i in range(len(mol_data.atomnum)):
        mol_data.atomnum[i] = i + start_val
    return mol_data

def renumber_residues(mol_data, start_val=1):
    '''
    Renumber residues so they start at 1
    '''
    reslist = []
    for i,resnum in enumerate(mol_data.resnum):
        if resnum not in reslist:
            newnum = len(reslist)
            reslist.append(resnum)
        else:
            newnum = reslist.index(resnum)
        mol_data.resnum[i] = newnum + start_val
    return mol_data

def rename_atoms(mol_data):
    '''
    Generate unique atom names
    '''
    for i,name in enumerate(mol_data.atomname):
    	mol_data.atomname[i] = ''.join([i.upper() for i in mol_data.atomname[i]
            if i.isalpha()])
    	
    occurrences = {}
    for i,atom in enumerate(mol_data.atomname):
    	if atom not in occurrences:
    		occurrences[atom] = [i,1]
    		mol_data.atomname[i] += str(occurrences[atom][1])
    	else:
    		occurrences[atom][1] += 1
    		mol_data.atomname[i] += str(occurrences[atom][1])
    return mol_data

def set_element(mol_data):
    '''
    Set atom element based on atom name
    '''
    for i,name in enumerate(mol_data.atomname):
        mol_data.element[i] = ''.join([i.title() for i in name if i.isalpha()])
    return mol_data

def has_hydrogen(mol_data):
    '''
    Check whether hydrogens are present
    '''
    return 'H' in mol_data.element.strip()

def strip_hydrogen(mol_data):
    '''
    Strip out all the hydrogens
    '''
    h_indices = [i for i,elem in enumerate(mol_data.element) if elem.strip() ==
            'H']
    new_mol_data = {}
    for key in mol_data.keys():
        new_mol_data[key] = [mol_data[key][i] for i in
                range(len(mol_data[key])) if i not in h_indices]
    return new_mol_data

def is_protein(mol_data):
    aa = standard_res.intersection(mol_data.resname)
    return aa > 0

def writepdb(mol_list, fname):
    '''
    Write molecule data to a file
    '''
    with open(fname, 'w') as f:
        chains_seen = []
        for mol_data in mol_list:
            mol_data_list = []
            for i in range(len(mol_data.atomnum)):
                j = 0
                for fieldwidth in pdb_fieldwidths:
                    if fieldwidth > 0:
                        fieldname = pdb_fieldnames[j]
                        mol_data_list[i].append('{:>{}s}'.format(str(mol_data[fieldname][i]),fieldwidth))
                        j += 1
                    else:
                        mol_data_list[i].append('{:>{}s}'.format('',abs(fieldwidth)))
                if mol_data['chainid'] not in chains_seen:
                    f.write('TER\n')
                    chains_seen.append(mol_data['chainid'])
                f.write(''.join(mol_data[i]))
                f.write('\n')
            f.write('TER\n')
        f.write('END\n')
