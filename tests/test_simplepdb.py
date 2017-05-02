import unittest
import env
import simplepdb as pdb

class IOTests(unittest.TestCase):
    '''
    Tests PDB file parsing and writing.
    '''
    def setUp(self):
        self.pdb = pdb.parse_pdb('3EML.pdb','r')

    def test_natoms(self):
        self.failUnless(1)

    def test_ters(self):
        self.failUnless(1)

    def test_data(self):
        self.failUnless(1)

class ActionTests(unittest.TestCase):
    '''
    Tests methods that manipulate molecular data.
    '''
    def setUp(sefl):
        self.receptor = pdb.parse_pdb('receptor.pdb', 'r')
        self.ligand = pdb.parse_pdb('ZMA.pdb', 'r')

    def test_get_ters(self):
        self.failUnless(1)

if __name__ == '__main__':
    unittest.main()
