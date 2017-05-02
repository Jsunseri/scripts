#!/usr/bin/python
import unittest
import env
import simplepdb as pdb
import os
from plumbum.cmd import awk, head, tail

class IOTests(unittest.TestCase):
    '''
    Tests PDB file parsing and writing.
    '''
    def setUp(self):
        self.pdb = pdb.simplepdb('chignolin.pdb')
        self.maxDiff = None

    def test_natoms(self):
        self.assertEqual(self.pdb.natoms, 138)

    def test_ters(self):
        self.assertEqual(self.pdb.ters[0], 10)

    def test_recordname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,0,6)}'])()
        #plumbum converts encoding to unicode and also adds an empty line at
        #the end...
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['recordname'], out[:-1])

    def test_atomnum(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,7,5)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['atomnum'], out[:-1])

    def test_atomname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,13,4)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['atomname'], out[:-1])

    def test_altloc(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,17,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['altloc'], out[:-1])

    def test_resname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,18,3)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['resname'], out[:-1])

    def test_chainid(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,22,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['chainid'], out[:-1])

    def test_resnum(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,23,4)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['resnum'], out[:-1])

    def test_rescode(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,27,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['rescode'], out[:-1])

    def test_x(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,31,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['x'], out[:-1])

    def test_y(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,39,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['y'], out[:-1])

    def test_z(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,47,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['z'], out[:-1])

    def test_occupancy(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,55,6)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['occupancy'], out[:-1])

    def test_beta(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,61,6)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['beta'], out[:-1])

    def test_element(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,77,2)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['element'], out[:-1])

    def test_charge(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,79,2)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['charge'], out[:-1])

    def test_write(self):
        fname = 'tmp.pdb'
        self.pdb.writepdb(fname)
        new_pdb = pdb.simplepdb(fname)
        os.remove(fname)
        self.assertEqual(self.pdb, new_pdb)

class ActionTests(unittest.TestCase):
    '''
    Tests methods that manipulate molecular data.
    '''
    def setUp(self):
        self.receptor = pdb.simplepdb('receptor.pdb')
        self.ligand = pdb.simplepdb('ZMA.pdb')

    def test_get_ters(self):
        self.failUnless(1)

if __name__ == '__main__':
    unittest.main()
