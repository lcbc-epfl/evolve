#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

Unittests for the gaapi modules

Tests
- crossovers
- generators
- Inidividual
- mutators
- replacers
- selectors

'''

import unittest
import openbabel
#from openbabel import openbabel
from  src import JobConfig
from src.gaapi import operator_types
from src.gaapi import generators
from src.gaapi import Individual
import main as main_program
import shutil
import os


def setUpModule():
    pass


def tearDownModule():
    pass



class TestGAAPI(unittest.TestCase):
    """Unittest."""

    def setUp(self):
        """Method to prepare the test fixture. Run BEFORE the test methods."""

        self.settings = JobConfig.Settings('example_files/test_individual.in')
        ga = main_program.initialise_ga(self.settings)
        self.individualone= Individual.Individual(self.settings)
        self.individualtwo = Individual.Individual(self.settings)
        pass

    def tearDown(self):
        """Method to tear down the test fixture. Run AFTER the test methods."""
        pass

    def addCleanup(self, function, *args, **kwargs):
        """Function called AFTER tearDown() to clean resources used on test."""
        pass


    #@unittest.skip("do not execute")
    def test_comparison_onedim(self):
        #test dominiation for one fitness
        self.individualone.fitnesses=[-10]
        self.individualtwo.fitnesses=[10]
        self.assertLess(self.individualone, self.individualtwo)
        self.assertGreater(self.individualtwo, self.individualone)

    ##@unittest.expectedFailure  #multidim not implented yet
    #@unittest.skip("do not execute")
    def test_comparison_multidim(self):
        # test dominiation for multiple fitnesses
        self.individualone.fitnesses = [-10, -10, -10]
        self.individualtwo.fitnesses = [10, 10, 10]
        self.assertLess(self.individualone, self.individualtwo)
        self.assertGreater(self.individualtwo, self.individualone)

    def test_applyComposition(self):
        print('applying compostion')
        #print(composition_residue_indexes)
        print(self.settings)
        self.composition_residue_indexes=[2, 4]
        self.composition_lower_bounds=[4, 5]
        self.composition_upper_bounds=[5, 6]
        print(self.composition_residue_indexes)
        self.individualone.applyComposition(self.settings)
        self.individualone.saveMol('output_files/applied_composition.pdb')



if __name__.__contains__("__main__"):
    print(__doc__)
    unittest.main()
