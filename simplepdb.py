#!/usr/bin/python
import pdb_util as util

class simplepdb:
    '''
    Parses and writes PDB files, and exposes limited functionality for
    manipulating their contents with a particular focus on the kinds of
    manipulations required for setting up MD simulations. 

    Attributes:
        mol_data: Dictionary of PDB column names and their values.
        ters_recs : Locations of breaks in the molecule, per the input PDB or
        resulting from simple operations such as merging molecules. Specified
        as a list of residues that appear immediately _before_ a break.
        natoms: Number of atoms in molecule(s).
    '''
    def __init__(self, pdb):
        '''
        Return a simplepdb object created by parsing the input PDB file.
        Can't construct an object without an input file because no utilities are 
        provided that could be used to construct a reasonable molecule.
        '''
        self.mol_data = self.parse_pdb(pdb)
        self.ters = self.get_ters(pdb)
        self.natoms = len(self.mol_data['atomnum'])

    def parse_pdb(self, pdb):
        '''
        Return a dictionary of PDB column names and their values for ATOM and
        HETATM records in the provided PDB file
        '''
        parse = util.make_parser(util.pdb_fieldwidths)
        f = open(pdb, 'r')
        mol_data_list = [parse(line) for line in f if line.startswith('HETATM') or line.startswith('ATOM')]
        f.close()
        mol_data = {}
        for i,field in enumerate(util.pdb_fieldnames):
            if (i == 1 or i == 6 or i == 8 or i == 9 or i == 10 or i == 11 or i ==
            12 or i == 14):
                fieldlist = [float(line[i]) for line in mol_data_list if
                        line[i].strip()]
            else:
                fieldlist = [line[i] for line in mol_data_list]
            mol_data[field] = fieldlist
        if not mol_data['element']:
            mol_data = set_element(mol_data)
        return mol_data

    def get_ters(self, pdb):
        '''
        Returns a list of breaks in a PDB file
        '''
        ters = []
        last_line = ''
        with open (pdb,'r') as f:
            for line in f:
                if line.startswith('TER'):
                    ter = float(line[22:26].strip())
                    if not ter:
                        ter = last_line[22:26].strip()
                    if ter: ters.append(ter) 
                last_line = line
        return ters

    def group_by_residue(self, mol_data):
        '''
        Rearrange atoms in a file so that atoms in the same residue are
        contiguous
        '''
        resmap = {}
        for old_idx in range(len(mol_data['atomnum'])):
            resnum = mol_data['resnum'][old_idx]
            if resnum not in resmap.keys():
                resmap[resnum] = [old_idx]
            else:
                resmap[resnum].append(old_idx)
        atnum = 1
        new_indices = [None] * len(mol_data['atomnum'])
        for res in resmap.keys():
            for old_idx in resmap[res]:
                mol_data['atomnum'][old_idx] = atnum
                new_indices.append(atnum-1)
                atnum += 1
        new_mol_data = {}
        for key in mol_data.keys():
            new_mol_data[key] = [mol_data[key][i] for i in new_indices]

        self.mol_data = new_mol_data
    
    def renumber_atoms(self, start_val=1):
        '''
        Renumber atoms so they start at 1
        '''
        for i in range(self.natoms):
            self.mol_data['atomnum'][i] = i + start_val
    
    def renumber_residues(self, start_val=1):
        '''
        Renumber residues so they start at 1
        '''
        reslist = []
        for i,resnum in enumerate(self.mol_data['resnum']):
            if resnum not in reslist:
                newnum = len(reslist)
                reslist.append(resnum)
            else:
                newnum = reslist.index(resnum)
            self.mol_data['resnum'][i] = newnum + start_val
    
    def rename_atoms(self):
        '''
        Generate unique atom names
        '''
        for i,name in enumerate(self.mol_data['atomname']):
        	self.mol_data['atomname'][i] = ''.join([i.upper() for i in
                    self.mol_data['atomname'][i]
                if i.isalpha()])
        	
        occurrences = {}
        for i,atom in enumerate(self.mol_data['atomname']):
        	if atom not in occurrences:
        		occurrences[atom] = [i,1]
        		self.mol_data['atomname'][i] += str(occurrences[atom][1])
        	else:
        		occurrences[atom][1] += 1
        		self.mol_data['atomname'][i] += str(occurrences[atom][1])
    
    def set_element(self):
        '''
        Set atom element based on atom name
        '''
        for i,name in enumerate(self.mol_data['atomname']):
            self.mol_data['element'][i] = ''.join([i.title() for i in name if i.isalpha()])

    def sanitize(self):
        '''
        Perform atom renumbering, residue renumbering, and regrouping atoms so
        residues are contiguous; if a small molecule, also uniquely names atoms
        and sets the element field
        '''
        group_by_residue()
        renumber_atoms()
        renumber_residues()
        if not is_protein():
            rename_atoms()
            set_element()

    def writepdb(self, fname, end=True):
        '''
        Write molecule data to a file
        '''
        self.sanitize()
        with open(fname, 'a') as f:
            for i in range(self.natoms):
                for fieldwidth in util.pdb_fieldwidths:
                    if fieldwidth > 0:
                        fieldname = util.pdb_fieldnames[j]
                        f.write('{:>{}s}'.format(str(self.mol_data[fieldname][i]),fieldwidth))
                    else:
                        f.write('{:>{}s}'.format('',abs(fieldwidth)))
                f.write('\n')
                if (self.mol_data['resnum'][i] in self.ters and
                self.mol_data['resnum'][i] != self.mol_data['resnum'][i+1]):
                    f.write('TER\n')
            if end==True: f.write('END\n')

    def has_hydrogen(self):
        '''
        Returns true if hydrogens are present
        '''
        return 'H' in self.mol_data['element'].strip()
    
    def strip_hydrogen(self):
        '''
        Strip out all the hydrogens
        '''
        h_indices = [i for i,elem in enumerate(self.mol_data['element']) if elem.strip() ==
                'H']
        new_mol_data = {}
        for key in self.mol_data.keys():
            new_mol_data[key] = [self.mol_data[key][i] for i in
                    range(len(self.mol_data[key])) if i not in h_indices]
        self.mol_data = new_mol_data
    
    def is_protein(self, ff=''):
        '''
        Returns true if standard amino acid residues are present
        '''
        aa = util.get_available_res(ff).intersection(self.mol_data['resname'])
        return len(aa) > 0
    
