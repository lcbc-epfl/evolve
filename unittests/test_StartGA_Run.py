#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

Unittests for starting a job

'''

import unittest
import openbabel
from  src import JobConfig
from src.gaapi import operator_types
from src.gaapi import generators
from src.gaapi import Individual
import main as main_program
import shutil
import os
from pathlib import Path

def setUpModule():
    pass


def tearDownModule():
    try:
        os.remove("mdinfo")
        os.remove("logfile")
        os.remove("leap.log")
        os.remove("new_mol.pdb")
        shutil.rmtree("amber_run")
    except:
        print("One of temporary files written by AMBER could not be deleted.")
    pass



class TestJobConfig(unittest.TestCase):
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

    def test_general_loading(self):
        print('Running tests for general parsing of input file')
        print("-" * 50)
        settings = JobConfig.Settings('example_files/ga_input.in')

        self.assertEqual(settings.verbose, True)
        self.assertEqual(settings.composition_optimization, True)
        self.assertEqual(len(settings.composition_lower_bounds), 2)
        self.assertEqual(len(settings.composition_lower_bounds), len(settings.composition_upper_bounds))
        self.assertEqual(settings.selector, operator_types.TOURNAMENT_WOR )
        self.assertEqual(settings.replacer, operator_types.ELITIST)

    def test_initPop_helicalstability(self):
        print('Running tests for initalizing of population for evaluator helical stability')
        print("-"*50)
        settings = JobConfig.Settings('example_files/ga_input.in')
        ga =main_program.initialise_ga(settings)
        parent_population= generators.initialisePopulation(settings)
        self.assertEqual(len(parent_population), 10)
        self.assertIsInstance(parent_population[0], Individual.Individual )

        generators.unbiased_protein_composition_generator(settings, parent_population)

        for i, indiv in enumerate(parent_population):
            self.assertEqual(indiv.composition[1], 175)
            self.assertEqual(len(indiv.fitnesses), len(settings.evaluators))
            self.assertEqual(len(indiv.composition), 2)
            # TODO
            # Test for length of indiv.mol
        
        # do an example fitness calculation of a defined individual = initial fitness test
        initial_indiv = Individual.Individual(settings)
        settings.originalResidues=['ALA', 'ALA']
        ga["evaluators"][0](settings, [initial_indiv], 0)
        self.assertAlmostEqual(initial_indiv.fitnesses[0], -64.173)
        f = open("amber_run/amber.out", 'r')
        for lineno, line in enumerate(f):
            if lineno==1604:
                final_energy=float(line[10:26])
        f.close()
        
        self.assertEqual(initial_indiv.fitnesses[0], final_energy)
        # check that new molecule file is identical to minimized file
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("pdb")
        obConversion.WriteFile(initial_indiv.mol, "new_mol.pdb")

        def read_atom_lines(infile):
            out=''
            for lineno, line in enumerate(infile):
                if line[:4] == "ATOM":
                    out+=line
            return out

        f1=open("new_mol.pdb", 'r')
        f2=open('amber_run/min_struct.pdb', 'r')
        f1content=read_atom_lines(f1)
        f2content = read_atom_lines(f2)
        self.assertTrue(f1content==f2content)
        f1.close()
        f2.close()


        # set predefined fitnesses and compare individuals to compare whether replacing works

    def test_initPop_helicalstability_openmm(self):
        print('OpenMM test')
        print("-"*50)
        settings = JobConfig.Settings('example_files/openmm.in')
        ga=main_program.initialise_ga(settings)
        parent_population= generators.initialisePopulation(settings)
        self.assertEqual(len(parent_population), 10)
        self.assertIsInstance(parent_population[0], Individual.Individual )

        generators.unbiased_protein_composition_generator(settings, parent_population)

        for i, indiv in enumerate(parent_population):
            self.assertEqual(indiv.composition[1], 175)
            self.assertEqual(len(indiv.fitnesses), len(settings.evaluators))
            self.assertEqual(len(indiv.composition), 2)
         
        # do an example fitness calculation of a defined individual = initial fitness test
        initial_indiv = Individual.Individual(settings)
        settings.originalResidues=['ALA', 'ALA'] # need to set this here, otherwise done in evaluator method.
        settings.gpu_openmm=True

        ga["evaluators"][0](settings, [initial_indiv], 0)
        print("openMM", initial_indiv.fitnesses[0])



if __name__.__contains__("__main__"):
    print(__doc__)
    unittest.main()
