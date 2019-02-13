'''
Created on Oct 5, 2016

@author: nbrownin
'''
import numpy as np
import ConfigParser 
from gaapi import operator_types 
import sys
import openbabel
import collections


class Settings(object):
    '''
    classdocs
    '''
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
        else:
            self.composition_optimization = False
        
        if (self.composition_optimization):
            self.composition_residue_indexes = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_residue_indexes').split()]
            self.composition_lower_bounds = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_lower_bounds').split()]
            self.composition_upper_bounds = [int(v) for v in self.config.get('OPTIMIZATION', 'composition_upper_bounds').split()]
            self.composition_library = self.config.get('OPTIMIZATION', 'composition_library')
        
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
            if (self.mutator == operator_types.POLY):
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
            
        # PARSE EVALUATORS
        if (self.config.has_section('EVALUATOR')):
            self.evaluators = [str(v) for v in self.config.get('EVALUATOR', 'evaluators').split()]
            if "turbomole_scf_energy" in self.evaluators:
                # parse turbomole specifics
                self.turbomole_template = self.config.get('EVALUATOR', 'turbomole_template')
            if "amber" in self.evaluators:
                # parse amber computation parameters
                # self.amber_path = self.config.get('EVALUATOR', 'amber_path')
                # self.mpi_path = self.config.get('EVALUATOR', 'mpi_path')
                self.use_gaussian = self.config.getboolean('EVALUATOR', 'use_gaussian')
                self.tleap_template = self.config.get('EVALUATOR', 'tleap_template')
                self.amber_params = self.config.get('EVALUATOR', 'amber_params')
                self.mpi_procs = int(self.config.get('EVALUATOR', 'mpi_processors'))
    
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
            
