'''

optimise.py

@author: Nicholas Browning
'''

import numpy as np
import matplotlib.pyplot as plt
from src.montecarlo import dihedral_opt as mcopt

'''
general, gly, prepro, pro
'''

import os
import sys

def read_input_file(inputFilePath):
    import ConfigParser
    settings = ConfigParser.ConfigParser()
    settings.read(inputFilePath)
    return settings

def doMain():
    import argparse
    from distutils.util import strtobool

    parser = argparse.ArgumentParser()
    parser.add_argument("-input_file", help="input file for montecarlo dihedral optimisation", type=str)
    args = parser.parse_args()
    
    settings = read_input_file(args.input_file)
    print  settings
    residue_pointers = [int (v) for v in  settings.get('SETTINGS', 'RESIDUE_POINTERS').split()]
    plot_optimised_distribution = strtobool(settings.get('SETTINGS', 'PLOT_OPTIMISED_DISTRIBUTION'))
    
    distribution_path = settings.get('DISTRIBUTION', 'DISTRIBUTION_FILES')
        
    probs = mcopt.loadDefaultProbabilityDistributions(distribution_paths=distribution_path)
    
    num_phipsi = len(residue_pointers)
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi, 2))
    
    print "Initial Solution:"
    print initial_solution
    
    move_lows = [-180] * num_phipsi
    move_highs = [180] * num_phipsi
    
    final_solution = mcopt.optimisePsiPhi(probs, residue_pointers, initial_solution, int(settings.get('SETTINGS', 'NUM_MC_STEPS')), move_lows, move_highs)
    
    print "Final Solution:"
    print final_solution
    
    
    if (plot_optimised_distribution):
        mcopt.plotProbabilityDistributionAndSolution(probs, residue_pointers, final_solution)
    
    if settings.has_section('MOLBUILD'):
        from src.gaapi import RotamerCreator as rc
        molbuild_residue_numbers = settings.get('MOLBUILD', 'MOLBUILD_RESIDUE_NUMBERS').split()
        molbuild_input = settings.get('MOLBUILD', 'MOLBUILD_INPUT')
        molbuild_output = settings.get('MOLBUILD', 'MOLBUILD_OUTPUT')
        # print molbuild_residue_numbers
        dihedral_data = np.zeros(len(final_solution), dtype='int32, float32, float32')
        # print dihedral_data
        for i in xrange (0, len(dihedral_data)):
            dihedral_data[i][0] = int(molbuild_residue_numbers[i])
            dihedral_data[i][1] = final_solution[i][0]
            dihedral_data[i][2] = final_solution[i][1]
        
        rc.buildMolBuildInput(input_structure_path=molbuild_input, output_structure_path=molbuild_output, output_path='molbuild.in', mode=rc.MOL_BUILD_MODE_DIHEDRAL_ONLY, dihedral_data=dihedral_data)
        print "Created MoleculeBuilder Input: ", 'molbuild.in'
if __name__ == '__main__':
    doMain()
