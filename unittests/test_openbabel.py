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
#from openbabel import openbabel
from  src import MoleculeInfo as mi
from  src import MoleculeCreator as mc
import filecmp
import os
import subprocess

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
        try:
            os.remove("mdinfo")
            os.remove("logfile")
            os.remove("leap.log")
            os.remove("new_mol.pdb")
        except:
            print("One of temporary files written by AMBER could  not be deleted.")
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


    def test_MoleculeCreator_rotname(self):

        # TODO
        # rotatematrix
        class DemoSettings(object):
            pass
        settings=DemoSettings()
        settings.use_res_type = False
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "example_files/ALA20.pdb")

        frag = openbabel.OBMol()
        obConversion.SetInAndOutFormats("mol2", "pdb")
        obConversion.ReadFile(frag, "../share/HID/HD3.mol2")
        
        mc.swapsidechain(settings, mol, 1, frag)

        obConversion.WriteFile(mol, "output_files/test_swap.pdb")

        self.assertTrue(filecmp.cmp('output_files/test_swap.pdb', 'example_files/H1A19.pdb', shallow=False), 'MoleculeCreator.swapsidechain has not yielded the same file')
        try:
            os.remove("output_files/test_swap.pdb")
        except:
            print("test_swap could not be deleted.")

    def test_MoleculeCreator_resname(self):

        # TODO
        # rotatematrix
        class DemoSettings(object):
            pass
        settings=DemoSettings()
        settings.use_res_type = True
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "example_files/ALA20.pdb")
 
        frag = openbabel.OBMol()
        obConversion.SetInAndOutFormats("mol2", "pdb")
        obConversion.ReadFile(frag, "../share/HID/HD3.mol2")

        mc.swapsidechain(settings, mol, 1, frag)

        obConversion.WriteFile(mol, "output_files/test_swap_resname.pdb")

        self.assertTrue(filecmp.cmp('output_files/test_swap_resname.pdb', 'example_files/H1A19_resname.pdb', shallow=False), 'MoleculeCreator.swapsidechain has not yielded the same file')
        try:
            os.remove("output_files/test_swap_resname.pdb")
        except:
            print("test_swap could not be deleted.")
    
    @unittest.skip("see if this one fails")
    def test_MoleculeCreator_glycine(self):

        class DemoSettings(object):
            pass
        settings=DemoSettings()
        settings.use_res_type = False
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol2", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "example_files/test_gly.mol2")

        frag = openbabel.OBMol()
        obConversion.SetInAndOutFormats("mol2", "pdb")
        obConversion.ReadFile(frag, "../share/ALA/A01.mol2")
        print('test')
        mc.swapsidechain(settings, mol, 8, frag)

        obConversion.WriteFile(mol, "output_files/test_gly_to_K.pdb")

        
        for rot_or_type in [False, True]:
            settings.use_res_type = rot_or_type
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, "example_files/3GB1WT_2clust2.mol2")
            frag = openbabel.OBMol()
            obConversion.ReadFile(frag, "../share/GLY/G01.mol2")
            mc.swapsidechain(settings, mol, 1, frag)
            obConversion.WriteFile(mol, "output_files/test_thr_to_gly_restype"+str(rot_or_type)+".pdb")
            obConversion.WriteFile(mol, "output_files/mol.pdb")
            os.chdir("./output_files/")
            print(os.getcwd())
            subprocess.call("tleap -f ../example_files/leap_glycine.in", shell=True)
            WasInputFileGenerated=os.path.isfile("mol.inpcrd")
            if WasInputFileGenerated:
                os.remove("mol.inpcrd")
            os.chdir("../")
            print("For ResType =", rot_or_type, "the result was ",WasInputFileGenerated)
            self.assertTrue(WasInputFileGenerated)

    def test_MoleculeCreator_histidine(self):

        class DemoSettings(object):
            pass
        settings=DemoSettings()
        settings.use_res_type = False
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol2", "pdb")

        for his in [("HID", "HD1"),("HID", "HD2"),("HID", "HD3"), ("HID", "HD6"), ("HID", "HD5"), ("HID", "HD8"), ("HIE", "HE2"), ("HIE", "HE7"), ("HIP", "HP7"), ("HIP", "HP3")]:
            for rot_or_type in [False, True]:
                settings.use_res_type = rot_or_type
                mol = openbabel.OBMol()
                obConversion.ReadFile(mol, "example_files/3GB1WT_2clust2.mol2")
                frag = openbabel.OBMol()
                obConversion.ReadFile(frag, "../share/"+his[0]+"/"+his[1]+".mol2")
                mc.swapsidechain(settings, mol, 1, frag)
                obConversion.WriteFile(mol, "output_files/test_thr_to_his_restype"+str(rot_or_type)+".pdb")
                obConversion.WriteFile(mol, "output_files/mol.pdb")
                os.chdir("./output_files/")
                print(os.getcwd())
                subprocess.call("tleap -f ../example_files/leap_glycine.in", shell=True)
                WasInputFileGenerated=os.path.isfile("mol.inpcrd")
                if WasInputFileGenerated:
                    os.remove("mol.inpcrd")
                os.chdir("../")
                print("For "+his[0]+" "+his[1]+" ResType =", rot_or_type, "the result was ",WasInputFileGenerated)
                self.assertTrue(WasInputFileGenerated)
        # self.assertTrue(filecmp.cmp('test_swap.pdb', 'example_files/H1A19.pdb', shallow=False), 'MoleculeCreator.swapsidechain has not yielded the same file')
        # try:
        #     os.remove("test_swap.pdb")
        # except:
        #     print("test_swap could not be deleted.")



if __name__.__contains__("__main__"):
    print(__doc__)
    unittest.main()
