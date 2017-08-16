'''
main.py

@author: Nicholas Browning
'''

from src.gaapi import selectors
from src.gaapi import mutators
from src.gaapi import crossovers
from src.gaapi import generators
from src.gaapi import replacers
from src import JobConfig
import evaluators
import argparse

import numpy as np

import printResidueInfo as pri

import openbabel

chosenSelector = selectors.tournamentSelectionWOR
chosenMutator = mutators.polynomial_dihedral
chosenCrossover = crossovers.simulated_binary_dihedral
chosenReplacer = replacers.elitist_replace
evaluators = [evaluators.testEnergyByMinimisation]


parent_population = []

def mainLoop(settings):
    curr_iteration = -1
    
    if settings.seed_population:
        # TODO - IMPLEMENT
        pass
    else:
        # generateInitialPopulation
    
        parent_population = generators.initialisePopulation(settings)
    
        if (settings.composition_optimization):
            print "Generating initial population composition from uniform distribution"
            generators.UniformCompositionGenerator(settings, parent_population)
            pass
        
        if (settings.dihedral_optimization):
            if (settings.mc_generate_dihedrals):
                print "Generating Initial Population dihedrals via MC"
                generators.MonteCarloDihedralGenerator(settings, parent_population)
            else:
                print "Generating Initial Population dihedrals from uniform distribution"
                generators.uniformDihedralGenerator(settings, parent_population)
        
        for i, eval in enumerate(evaluators):
            eval(settings, parent_population, i)
        
    print "Iterating"
    curr_iteration = 0
    
    while (curr_iteration < settings.max_iteration):
        
        
        printIterationInfo(settings, curr_iteration, parent_population)
        # select
        indexes = selectors.selector(settings, parent_population, chosenSelector)
        
        if (settings.verbose):
            print "Mating pool:", indexes
            
        child_population = [parent_population[i] for i in indexes]
        
        # crossover       
        child_population = crossovers.crossover(settings, child_population, chosenCrossover)
        
        # mutate
        child_population = mutators.mutate(settings, child_population, chosenMutator)
        
        
        # evaluate - build new mol files, run amber minimisation/whatever else, update individuals with resulting dihedrals
        for i, eval in enumerate(evaluators):
            eval(settings, child_population, i)
        # elitism
        
            
        parent_population = replacers.replace(settings, parent_population, child_population, chosenReplacer)
        
        curr_iteration += 1
    

    
def printIterationInfo(settings, curr_iteration, pop):
    
    print "Current Iteration: ", curr_iteration
    
    min_indiv = np.argmin([p.fitnesses[0] for p in pop])
    
    print "Population Average Fitness(es):", np.mean([p.fitnesses for p in pop])
    print "Population STD(s):", np.std([p.fitnesses for p in pop])
    print "best individual:", min_indiv
    
    if (pop[min_indiv].dihedral_residue_indexes != None):
        print pop[min_indiv].phi_dihedrals, pop[min_indiv].psi_dihedrals
        
    print "fitness(es):", pop[min_indiv].fitnesses
    
    
    if (settings.verbose):
        for i in xrange (settings.population_size):
            print i, pop[i], pop[i].fitnesses, pop[i].phi_dihedrals, pop[i].psi_dihedrals, pop[i].mol
            
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdb")

    obConversion.WriteFile(pop[min_indiv].mol, "min_" + str(curr_iteration) + ".pdb") 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-input_file', type=str)
    
    args = parser.parse_args()
    
    settings = JobConfig.Settings(args.input_file)
    print "-- START SETTINGS PRINT --"
    settings.printSettings()
    print "-- END SETTINGS PRINT -- "
    
    print "--- START MOLECULE PRINT --"
    pri.printResidueInfo(settings.initial_molecule)
    print "--- END MOLECULE PRINT --"
    
    print "Running GA"
    mainLoop(settings)
    print "Finished GA"
