'''
main.py

@author: Nicholas Browning
'''

import numpy as np
import argparse
import openbabel
import time
import json

from src import JobConfig
from src.gaapi import Individual
from src.gaapi import selectors
from src.gaapi import mutators
from src.gaapi import crossovers
from src.gaapi import generators
from src.gaapi import replacers
import evaluators as evals
import printResidueInfo as pri


def mainLoop(settings):
    
    curr_iteration = -1
    parent_population = []

    ga = initialise_ga(settings)

    # Compute fitness of initial structure given to GA and store in file
    print('Fitness of {}:'.format(settings.initial_molecule_path))
    initial_indiv = Individual.Individual(settings)
    initial_indiv.init = True
    
    for k, eval in enumerate(ga["evaluators"]):
        print("EVALUATORS {}".format(ga["evaluators"][k].__name__))
        eval(settings, [initial_indiv], k)

    with open(settings.output_path + '/initial_fitness.dat', 'w') as outfile:  
        json.dump(initial_indiv.fitnesses[0], outfile)

    parent_population = generators.initialisePopulation(settings)
    
    if (settings.dihedral_optimization):
        # generators.BasiliskSideChainDihedralsGenerator(settings, parent_population)
        generators.UniformDihedralGenerator(settings, parent_population)
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("pdb")
        for i in range(len(parent_population)):
            obConversion.WriteFile(parent_population[i].mol, settings.output_path + '/' + "initial_mol_" + str(i) + ".pdb")


    for i, eval in enumerate(ga["evaluators"]):
        print("EVALUATORS {}".format(ga["evaluators"][i].__name__))
        eval(settings, parent_population, i)
    
    
    
    curr_iteration = 0
    best_fitnesses = [np.min([p.fitnesses[0] for p in parent_population])]
 
    while (curr_iteration < settings.max_iteration):

        printIterationInfo(settings, curr_iteration, parent_population, False)
        
        print "\n------------Iterating"
        # select
        indexes = selectors.selector(settings, parent_population, ga["chosenSelector"])
        if (settings.verbose):
            print "Mating pool:", indexes  
        child_population = [parent_population[i] for i in indexes]
        
        # crossover
        child_population = crossovers.crossover(settings, child_population, ga["chosenCrossover"])
        
        # mutate
        child_population = mutators.mutate(settings, child_population, ga["chosenMutator"])
        
        # evaluate - build new mol files, run amber minimisation/whatever else, update individuals with resulting dihedrals
        for i, eval in enumerate(ga["evaluators"]):
            print("EVALUATORS {}".format(ga["evaluators"][i].__name__))
            eval(settings, child_population, i)

        # elitism 
        parent_population = replacers.replace(settings, parent_population, child_population, ga["chosenReplacer"])

        best_fitnesses.append(min([p.fitnesses[0] for p in parent_population]))
        curr_conv_relative_error = abs((best_fitnesses[curr_iteration + 1] - best_fitnesses[curr_iteration]) / best_fitnesses[curr_iteration + 1])
        print("Convergence relative error: {}".format(curr_conv_relative_error))

        curr_iteration += 1

        with open(settings.output_path + '/fitnesses.dat', 'w') as outfile:  
            json.dump(best_fitnesses, outfile)
        
    printIterationInfo(settings, curr_iteration, parent_population, True)
    print('Best GA fitnesses written in fitnesses.dat')

    
def initialise_ga(settings):
    '''
    chosenSelector = selectors.tournamentSelectionWOR
    chosenMutator = mutators.polynomial_dihedral
    chosenCrossover = crossovers.simulated_binary_dihedral
    chosenReplacer = replacers.elitist_replace
    evaluators = [evaluators.amber_energy]#[evaluators.testEnergyByMinimisation]
    '''

    chosenSelector, chosenMutator, chosenCrossover, chosenReplacer, evaluators = None, None, None, None, None
    
    if (settings.selector == "tournament_wor"):
        chosenSelector = selectors.tournamentSelectionWOR

    if (settings.mutator == "poly_mutation"):
        chosenMutator = mutators.polynomial_mutation
    elif (settings.mutator == "uniform_mutation"):
        chosenMutator = mutators.uniform_mutation
        
    if (settings.crossover == "sbx"):
        chosenCrossover = crossovers.simulated_binary_crossover
    elif (settings.crossover == "uniform_crossover"):
        chosenCrossover = crossovers.uniform_crossover
    
    if (settings.replacer == "elitism"):
        chosenReplacer = replacers.elitist_replace
    
    
    evaluators = []
    for k in settings.evaluators:
        if ("openbabel_ff94" in settings.evaluators):
            evaluators.append(evals.testEnergyByMinimisation)
        elif ("amber" in settings.evaluators):
            evaluators.append(evals.amber_energy_simplified)
        else:
            pass  # could add other choices in evaluators.py

    return {'chosenSelector': chosenSelector, 'chosenMutator': chosenMutator, 'chosenCrossover': chosenCrossover, 'chosenReplacer': chosenReplacer, 'evaluators': evaluators}



    
def printIterationInfo(settings, curr_iteration, pop, ending):

    print("CURRENT POPULATION")
    if (settings.basilisk_and_sidechains and settings.verbose):
        
        print("ind, ind fitness(es), res, phi, psi, chi angles")
        for i in xrange (settings.population_size):
            for j in xrange(0, len(pop[i].dihedral_residue_indexes)):
                print("{0}, {5}, {1}, {2}, {3}, {4}".format(i, pop[i].dihedral_residue_indexes[j], pop[i].phi_dihedrals[j], pop[i].psi_dihedrals[j], pop[i].chi_angles[j], pop[i].fitnesses))

    if (settings.verbose and not settings.basilisk_and_sidechains):
        print("ind, ind fitness(es), res, phi, psi")
        for i in xrange (settings.population_size):
            for j in xrange(0, len(pop[i].dihedral_residue_indexes)):
                print("{0}, {4}, {1}, {2}, {3}".format(i, pop[i].dihedral_residue_indexes[j], pop[i].phi_dihedrals[j], pop[i].psi_dihedrals[j], pop[i].fitnesses))

    print "Current Iteration: ", curr_iteration
    min_indiv = np.argmin([p.fitnesses[0] for p in pop])
    
    print "Population Average Fitness(es):", np.mean([p.fitnesses for p in pop])
    print "Population STD(s):", np.std([p.fitnesses for p in pop])
    print "Best individual:", min_indiv
    
    if (pop[min_indiv].dihedral_residue_indexes != None):
        print('Best ind - PHIs, PSIs')
        print pop[min_indiv].phi_dihedrals, pop[min_indiv].psi_dihedrals
        
    print "Best ind - Fitness(es):", pop[min_indiv].fitnesses
            
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdb")

    if (curr_iteration == settings.max_iteration or ending):
        print("OUT: FINAL MOLECULE in {}".format(settings.output_path + '/' + settings.output_file))
        obConversion.WriteFile(pop[min_indiv].mol, settings.output_path + '/' + settings.output_file)

    elif (settings.store_intermediate_mol):
        print("OUT: INTERMEDIATE MOLECULE in {}".format(settings.output_path + '/' + "min_mol_iter_" + str(curr_iteration) + ".pdb"))
        obConversion.WriteFile(pop[min_indiv].mol, settings.output_path + '/' + "min_mol_iter_" + str(curr_iteration) + ".pdb")



    
if __name__ == '__main__':
    start_time = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_file', type=str)
    args = parser.parse_args()

    settings = JobConfig.Settings(args.input_file)
    print("\n-- START SETTINGS PRINT --")
    settings.printSettings()
    print("-- END SETTINGS PRINT -- ")
    
    print("\n--- START MOLECULE PRINT --")
    pri.printResidueInfo(settings.initial_molecule)
    print("--- END MOLECULE PRINT --")
    
    print("\nRunning GA")
    mainLoop(settings)
    print("Finished GA in --- {:.2f} seconds ---".format(time.time() - start_time))
