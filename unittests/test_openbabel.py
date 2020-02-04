#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

Unittests for the openbabel based modules

Tests
- MoleculeCreator
- MoleculeInfo

'''

import unittest
import openbabel
from  src import MoleculeInfo as mi
from  src import MoleculeCreator as mc
import filecmp
import os

def setUpModule():
    pass


def tearDownModule():
    pass



class TestGAAPI(unittest.TestCase):
    """Unittest."""

    def setUp(self):
        """Method to prepare the test fixture. Run BEFORE the test methods."""
        pass

    def tearDown(self):
        """Method to tear down the test fixture. Run AFTER the test methods."""
        pass

    def addCleanup(self, function, *args, **kwargs):
        """Function called AFTER tearDown() to clean resources used on test."""
        pass

    @classmethod
    def setUpClass(cls):
        """Class method called BEFORE tests in an individual class run. """
        pass  # Probably you may not use this one. See setUp().

    @classmethod
    def tearDownClass(cls):
        """Class method called AFTER tests in an individual class run. """
        pass  # Probably you may not use this one. See tearDown().

    def test_MoleculeInfo(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "example_files/ALA20.pdb")
        obres = mol.GetResidue(0)
        self.assertTrue(obres.GetName()=="ACE")
        obres = mol.GetResidue(1)
        self.assertTrue(mi.getResType(obres) == 'ALA')
        #beta_atom = getBetaAtom(obres)
        # # TODO
        # test for proline and glycine

        # AlphaCarbon
        self.assertIs(mi.getAlphaCarbon(obres).GetIdx(),9)

        # bb carboxzl
        self.assertIs(mi.getBBCarboxyl(obres).GetIdx(),15)

        #bb nitrogen
        self.assertIs(mi.getBBNitrogen(obres).GetIdx(),7)


        # Dihedral tests

        # TODO
        self.assertAlmostEqual(mi.getChi1DihedralAngle(mol, obres),-58.54682580446737)

        x,y,z=mi.getPosition(mi.getAlphaCarbon(obres))
        self.assertTrue(x==-0.333)
        self.assertTrue(y==-0.155)
        self.assertTrue(z==-0.121)

        self.assertIs(mi.countBonds(mi.getAlphaCarbon(obres)),4)

    def test_MoleculeCreator(self):

        # TODO
        # rotatematrix

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "example_files/ALA20.pdb")

        frag = openbabel.OBMol()
        obConversion.SetInAndOutFormats("mol2", "pdb")
        obConversion.ReadFile(frag, "../share/HID/HD3.mol2")

        mc.swapsidechain(mol, 1, frag)

        obConversion.WriteFile(mol, "test_swap.pdb")

        self.assertTrue(filecmp.cmp('test_swap.pdb', 'example_files/H1A19.pdb', shallow=False), 'MoleculeCreator.swapsidechain has not yielded the same file')
        try:
            os.remove("test_swap.pdb")
        except:
            print("test_swap could not be deleted.")



if __name__.__contains__("__main__"):
    print(__doc__)
    unittest.main()
