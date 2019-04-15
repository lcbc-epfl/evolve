'''
Created on Oct 5, 2016

@author: nbrownin
'''
import argparse
import numpy as np
import ConfigParser 
from gaapi import operator_types 
import sys
import openbabel
import collections
from src import constants as cnts

class Settings(object):
    '''
    classdocs
    '''
    
    composition_optimization = False
    composition_residue_indexes = None
    composition_lower_bounds = None
    composition_upper_bounds = None
    composition_library = None
    allowed_residue_types = None
    
    
    backbone_dihedral_optimization = False
    dihedral_residue_indexes = None
    
    sidechain_dihedral_optimisation = False
    basilisk_and_sidechains = False

    solution_min_fitness = None

    seed_population = False
    population_file_path = None
    
    unbiased_protein_composition_generator = False
    
    population_size = 0
    max_iteration = 0
    initial_energy = 0.0
    
    helical_dielectric = 0
    
    initial_molecule = None
    
    def __init__(self, configFilePath):
        
        self.configFilePath = configFilePath
        
        self.config = ConfigParser.ConfigParser()
        self.config.read(configFilePath)
        
        # PARSE [MOLECULE]
        self.initial_molecule_path = str(self.config.get('MOLECULE', 'initial_molecule').strip('\''))
        # openbabel parse
        obConversion = openbabel.OBConversion()
        if (not obConversion.SetInFormat(self.initial_molecule_path.split(".")[1])):
            print "Problem reading ", self.initial_molecule_path, " filetype: ", self.initial_molecule_path.split(".")[1]
            sys.exit(0)
        print "Succesfully set molecule filetype:", self.initial_molecule_path.split(".")[1]
        
        self.initial_molecule = openbabel.OBMol()
        
        if (not obConversion.ReadFile(self.initial_molecule, self.initial_molecule_path)):
            print "Problem reading ", self.initial_molecule_path
            sys.exit(0)
        print "Succesfully read molecule:", self.initial_molecule_path
    
    
        # PARSE [OPTIMIZATION]
        if self.config.has_option('OPTIMIZATION', 'composition_optimization'):
            self.composition_optimization = self.config.getboolean('OPTIMIZATION', 'composition_optimization')

    
        if (self.composition_optimization):
            self.composition_residue_indexes = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_residue_indexes').split()]
            self.composition_lower_bounds = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_lower_bounds').split()]
            self.composition_upper_bounds = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_upper_bounds').split()]
            self.composition_library = self.config.get('OPTIMIZATION', 'composition_library')
           
            # Restrict energy/rotamer/types list to only those specified in input
            
        if (self.config.has_option('OPTIMIZATION', 'allowed_residue_types')):
            cnts.allowed_residue_types = [v.strip() for v in self.config.get('OPTIMIZATION', 'allowed_residue_types').split(' ')]
            cnts.selected_rotamers = cnts.subselect_rotamers(cnts.allowed_residue_types)
            cnts.selected_rotamer_types = cnts.subselect_rotamer_types(cnts.allowed_residue_types)
            print cnts.selected_rotamers, cnts.selected_rotamer_types
                
        if (self.config.has_option('OPTIMIZATION', 'min_fitness')):
            self.solution_min_fitness = data = float(self.config.get('OPTIMIZATION', 'min_fitness'))
        
            

        if self.config.has_option('OPTIMIZATION', 'backbone_dihedral_optimization'):
            self.backbone_dihedral_optimization = self.config.getboolean('OPTIMIZATION', 'backbone_dihedral_optimization')
        else:
            self.backbone_dihedral_optimization = False
    
        if (self.backbone_dihedral_optimization):
            self.dihedral_residue_indexes = [int(v) for v in self.config.get('OPTIMIZATION', 'dihedral_residue_indexes').split()]
            
            
        if self.config.has_option('OPTIMIZATION', 'sidechain_dihedral_optimisation'):
            self.sidechain_dihedral_optimisation = self.config.getboolean('OPTIMIZATION', 'sidechain_dihedral_optimisation')
        else:
            self.sidechain_dihedral_optimisation = False
        # PARSE [GA GLOBAL]
        if (self.config.has_option('GA_GLOBAL', 'seed_population')):
            self.seed_population = self.config.getboolean('GA_GLOBAL', 'seed_population')
            self.population_file_path = self.config.get('GA_GLOBAL', 'population_file_path')
        
        self.population_size = int(self.config.get('GA_GLOBAL', 'population_size'))
        self.max_iteration = int(self.config.get('GA_GLOBAL', 'max_iteration'))
        
        # PARSE [GA_SELECTION]
        if (self.config.has_section('GA_SELECTION')):
            self.selector = self.config.get('GA_SELECTION', 'selector')
            # print self.selector, operator_types.TOURNAMENT_WOR
            if (self.selector == operator_types.TOURNAMENT_WOR):
                self.tournament_size = int(self.config.get('GA_SELECTION', 'tournament_size'))
            
        # PARSE [GA_CROSSOVER]  
        if (self.config.has_section('GA_CROSSOVER')):
            self.crossover = self.config.get('GA_CROSSOVER', 'crossover')
            
            if (self.crossover == operator_types.SBX):
                self.sbx_distribution_index = float(self.config.get('GA_CROSSOVER', 'sbx_distribution_index'))
                
            self.mating_probability = float(self.config.get('GA_CROSSOVER', 'mating_probability'))
            self.genewise_crossover_probability = float(self.config.get('GA_CROSSOVER', 'genewise_crossover_probability'))
        
        # PARSE [GA MUTATION]
        if (self.config.has_section('GA_MUTATION')):
            self.mutator = self.config.get('GA_MUTATION', 'mutator')
            
            if (self.mutator == operator_types.POLY_M):
                self.poly_eta = float(self.config.get('GA_MUTATION', 'poly_eta'))
            
            self.mutation_probability = float(self.config.get('GA_MUTATION', 'mutation_probability'))
            self.genewise_mutation_probability = float(self.config.get('GA_MUTATION', 'genewise__mutation_probability'))
        
        # PARSE [GA GA_REPLACER]
        if (self.config.has_section('GA_REPLACER')):
            self.replacer = self.config.get('GA_REPLACER', 'replacer')
            if (self.replacer == operator_types.ELITIST):
                self.elitist_factor = float(self.config.get('GA_REPLACER', 'elitist_factor'))
            
        # PARSE [MC_GENERATE]
        if (self.config.has_section('MC_GENERATE')):
            self.mc_generate_dihedrals = self.config.getboolean('MC_GENERATE', 'mc_generate_dihedrals')
            self.dihedral_probability_pointers = [int(v) for v in self.config.get('MC_GENERATE', 'dihedral_probability_pointers').split()]
            self.distribution_path = self.config.get('MC_GENERATE', 'distribution_path')
            self.mcmove_lbound = float(self.config.get('MC_GENERATE', 'mcmove_lbound'))
            self.mcmove_ubound = float(self.config.get('MC_GENERATE', 'mcmove_ubound'))
            self.num_mc_steps = int(self.config.get('MC_GENERATE', 'num_mc_steps'))
    
        if (self.config.has_section('GENERATOR')):
            self.unbiased_protein_composition_generator = self.config.getboolean('GENERATOR', 'unbiased_protein_composition_generator')
        
        # PARSE EVALUATORS
        if (self.config.has_section('EVALUATOR')):
            self.evaluators = [str(v) for v in self.config.get('EVALUATOR', 'evaluators').split()]
            if "turbomole_scf_energy" in self.evaluators:
                # parse turbomole specifics
                self.turbomole_template = self.config.get('EVALUATOR', 'turbomole_template')
            if "amber" in self.evaluators:
                self.tleap_template = self.config.get('EVALUATOR', 'tleap_template')
                self.amber_params = self.config.get('EVALUATOR', 'amber_params')
                self.mpi_procs = int(self.config.get('EVALUATOR', 'mpi_processors'))
            if "helical_stability" in self.evaluators:
                self.tleap_template = self.config.get('EVALUATOR', 'tleap_template')
                self.amber_params = self.config.get('EVALUATOR', 'amber_params')
                self.mpi_procs = int(self.config.get('EVALUATOR', 'mpi_processors'))
                self.helical_dielectric = int(float(self.config.get('EVALUATOR', 'dielectric')))
                if (self.helical_dielectric == 80):
                    self.helical_dielectric = 2
                elif (self.helical_dielectric == 50):
                    self.helical_dielectric = 1
                elif (self.helical_dielectric == 10):
                    self.helical_dielectric = 0
                    
                
        # PARSE [IO]
        self.output_path = self.config.get('IO', 'output_path')
        self.output_file = self.config.get('IO', 'output_file')
        self.store_intermediate_mol = self.config.getboolean('IO', 'store_intermediate_mol')

        
        # PARSE [DEBUG]
        self.verbose = self.config.getboolean('DEBUG', 'verbose')
    

    def printSettings(self):
        ordered_dict = collections.OrderedDict(self.__dict__)  # No change since self.__dict__ already disordered..., could be removed
        for attribute, value in ordered_dict.items():
            print('{} : {}'.format(attribute, value))
            

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-input_file', type=str)
    args = parser.parse_args()
    settings = Settings(args.input_file)

