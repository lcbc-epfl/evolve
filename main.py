'''
This is the main module of EVOLVE.

It initialises the genetic algorithm search,
determines initial fitnesses for each evaluator and then loops over the generations.

.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu
'''
from __future__ import division
from __future__ import print_function

from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import argparse
import openbabel
import time
import json
import sys
import os
from src import constants as cnts
from src import JobConfig
from src.gaapi import Individual
from src.gaapi import selectors
from src.gaapi import mutators
from src.gaapi import crossovers
from src.gaapi import generators
from src.gaapi import replacers
from src import MoleculeInfo as mi
import src.outprocesses as op
from collections import Counter 
import evaluators as evals
# import  printResidueInfo as pri


def mainLoop(settings):
    '''

    The main GA loop that initalises the optimizaion procedure with an empty parent population.
    A seed population can be used (not implemented yet)
    Then the initial fitness is determined and a parent population is created (if not seeded)
    Its fitness is determined.

    In the main GA while loop, from each parent population a child generation is created
    through the crossover and mutation operations, their fitness is determined and
    the replacer operation is used to build the next parent generation.

    For single objective runs the comparison is made using simple sorting of the children and retention of a certain percentage of good solutions (elistism)
     For multiobjective runs Deb et.al. NSGA2 algorithm is used which autmatically incorporates elistim.

    If you want to add your own evaluation function you want to add your method to evaluators.py and add corresponding methods to outprocesses and methods.
    Your evaluator needs to be registered in operator_types.

    In order to bias the first parent population the generators file can be used.
    For composotion optimization we currently implement random aminoacid rotamer subsitution (favors aas with many rotamers), unform (equal probability for each allowed amino acid) and swissprot (based on swissprot probability densities)

    Parameters
    ----------
    settings: object
        see :class:`src.JobConfig.Settings`

    '''
    curr_iteration = -1
    parent_population = []

    ga = initialise_ga(settings)
    
    # if ("helical_stability" in settings.evaluators):
    #     print('Fitness of {}:'.format(settings.initial_molecule_path))
    #     initial_indiv = Individual.Individual(settings)
    #     initial_indiv.init = True
    #
    #     print("MINIMISING INITIAL STRUCTURE")
    #     for k, eval in enumerate(ga["evaluators"]):
    #         print("EVALUATORS {}".format(ga["evaluators"][k].__name__))
    #         eval(settings, [initial_indiv], k)
    #
    #     with open(settings.output_path + '/initial_fitness.dat', 'w') as outfile:
    #         json.dump(initial_indiv.fitnesses[0], outfile)
    #
    #     settings.initial_energy = initial_indiv.fitnesses[0]
    #
    #     print("INITIAL_ENERGY", settings.initial_energy)
    #
    #     directory = settings.output_path + "/amber_run"
    #
    #     obConversion = openbabel.OBConversion()
    #     obConversion.SetInFormat("pdb")
    #
    #     molec = openbabel.OBMol()
    #     obConversion.ReadFile(molec, directory + '/min_struct.pdb')
    #
    #     settings.initial_molecule = molec
    #
    if settings.seed_population:
        # TODO - IMPLEMENT
        pass
    else:
        # generateInitialPopulation
        if settings.initial_fitness_computation:
            # if multi individual we need to do this for every individual.
            if (settings.multi_individual):
                print('Fitness of initial frames:')
                pass
            else:
                print('Fitness of {}:'.format(settings.initial_molecule_path))
            if not os.path.exists(settings.output_path + "/work_dir"):
                os.makedirs(settings.output_path + "/work_dir")
                f = open(settings.output_path + "/work_dir/evaluator_mmpbsa_multi.log", "w+")
                f.close()
            initial_indiv = Individual.Individual(settings)
            initial_indiv.init = True
            if (settings.multi_individual):
                settings.originalResidues = [mi.getResType(initial_indiv.mol[0].GetResidue(j)) for j in settings.composition_residue_indexes]
            else:
                settings.originalResidues = [mi.getResType(initial_indiv.mol.GetResidue(j)) for j in settings.composition_residue_indexes]

            print("Original residues:", settings.originalResidues)
            print(ga["evaluators"])
            for k, eval in enumerate(ga["evaluators"]):
                print("EVALUATORS {}".format(ga["evaluators"][k].__name__))
                eval(settings, [initial_indiv], k)
                print("evaluated", ga["evaluators"][k])

            with open(settings.output_path + '/initial_fitness.dat', 'w') as outfile:
                json.dump(str(initial_indiv.fitnesses), outfile)

            if (settings.multi_individual):
                settings.initial_energy = initial_indiv.fitnesses[1]  # TODO only works for multi mmpbsa for now
                pass
            else:
                settings.initial_energy = initial_indiv.fitnesses[0]


            print("INITIAL_ENERGY", settings.initial_energy)

            pass

        parent_population = generators.initialisePopulation(settings)


        # TODO
        # This should become similar to the way evaluators are chosen
        if (settings.composition_optimization):
            generators.generate(settings, parent_population, ga["chosenGenerator"])
            #if (settings.unbiased_protein_composition_generator):
            #     print("Generating initial population composition from unbiased distribution")
            #     generators.unbiased_protein_composition_generator(settings, parent_population)
            #elif (settings.swissprot_composition_generator):
            #    print("Generating initial population composition from swissprot probability distribution")
            #    generators.swissprot_composition_generator(settings, parent_population)
            #else:
            #    print("Generating initial population composition from uniform distribution")
            #    generators.UniformCompositionGenerator(settings, parent_population)

        if (settings.multi_individual):
            print("Resetting logfile in work dir\n "+settings.output_path + "/work_dir/evaluator_mmpbsa_multi.log")
            if os.path.exists(settings.output_path + "/work_dir"):
                try:
                    os.remove(settings.output_path + "/work_dir"+"/evaluator_mmpbsa_multi.log")
                except OSError as e: # this would be "except OSError, e:" before Python 2.6
                    os.makedirs(settings.output_path + "/work_dir")
                    f= open(settings.output_path + "/work_dir/evaluator_mmpbsa_multi.log","w+")
                    f.close()
                    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
                        raise # re-raise exception if a different error occurred
                    pass

        if (settings.use_compute_cluster):
            print("Create folder for SLURM output\n "+settings.output_path + "/slurm")
            if not os.path.exists(settings.output_path +"/slurm"):
                os.makedirs(settings.output_path +"/slurm")

        if (settings.backbone_dihedral_optimization):
            if (settings.mc_generate_dihedrals):
                print("Generating Initial Population dihedrals via MC")
                generators.MonteCarloDihedralGenerator(settings, parent_population)

            else:
                print("Generating Initial Population dihedrals from uniform distribution")
                generators.UniformDihedralGenerator(settings, parent_population)
                
            if (settings.basilisk_and_sidechains):
                print("Backbone-dependent side chains generated by Basilisk")
                generators.BasiliskSideChainDihedralsGenerator(settings, parent_population)

        for i, eval in enumerate(ga["evaluators"]):
            print("EVALUATORS {}".format(ga["evaluators"][i].__name__))
            eval(settings, parent_population, i)
      
    curr_iteration = 0
    # If multiobjevtive optimization the population will already be sorted, the fittest individual is always 0
    # For single objective the fittest individual is not 0 and needs to be identified
    best_fitnesses = np.zeros(shape=(settings.max_iteration,1))
    if settings.no_evaluators==1:
        best_fitnesses = [np.min([p.fitnesses[0] for p in parent_population])]
    else:
        best_fitnesses=[parent_population[0].fitnesses]
        # TODO Need to implement convergence relative error for multi objective
 
    while (curr_iteration < settings.max_iteration):

        printIterationInfo(settings, curr_iteration, parent_population, False)
        
        print("\n------------Iterating")
        # select
        indexes = selectors.selector(settings, parent_population, ga["chosenSelector"])
        if (settings.verbose):
            print("Mating pool:", indexes)  
        child_population = [parent_population[i] for i in indexes]
        
        # crossover
        child_population = crossovers.crossover(settings, child_population, ga["chosenCrossover"])
        
        # mutate
        child_population = mutators.mutate(settings, child_population, ga["chosenMutator"])
        
        # evaluate - build new mol files, run amber minimisation/whatever else, update individuals with resulting dihedrals
       
        for j in range(0, len(child_population)):
            if (settings.multi_individual):
                child_population[j].mol = [None] * settings.no_frames
                for frame, x in enumerate(settings.initial_molecule):
                    child_population[j].mol[frame] = openbabel.OBMol(x)
                    pass
                pass
            else:
                child_population[j].mol = openbabel.OBMol(settings.initial_molecule)
                pass
            child_population[j].applyComposition(settings)
      
        for i, eval in enumerate(ga["evaluators"]):
            print("EVALUATORS {}".format(ga["evaluators"][i].__name__))
            eval(settings, child_population, i)

        # elitism 
        parent_population = replacers.replace(settings, parent_population, child_population, ga["chosenReplacer"])

        if settings.no_evaluators == 1:
            best_fitnesses.append(min([p.fitnesses[0] for p in parent_population]))
            curr_conv_relative_error = abs(
                (best_fitnesses[curr_iteration + 1] - best_fitnesses[curr_iteration]) / best_fitnesses[
                    curr_iteration + 1])
            print("Convergence relative error: {}".format(curr_conv_relative_error))
            pass
        else:
            #best_fitnesses.append(parent_population[0].fitnesses)
            #print("Best Individual fitness", parent_population[0].fitnesses)
            # TODO Need to implement convergence relative error for multi objective
            pass

        curr_iteration += 1
        settings.curr_iteration = curr_iteration
        
        if settings.no_evaluators == 1:
            np.savetxt(settings.output_path + '/fitnesses.dat', best_fitnesses, delimiter='\n')
            if list(best_fitnesses).count(best_fitnesses[-1])>settings.convergence_cycles:
                print(f"Reached convergence. Last {list(best_fitnesses).count(best_fitnesses[-1])} generations were identical.")
                printIterationInfo(settings, curr_iteration, parent_population, True)
                sys.exit(0)
        else:
            # in order to correctly write out with np savetext we need to reduce multidimensional fitness arrays to a 1d list
            fitness_list = []
            fitness_line = ''
            for fitness in best_fitnesses:
                for fitness_ind in fitness:
                    fitness_line += str(fitness_ind) + " "
                pass
                fitness_list.append([fitness_line])
                fitness_line = ''
            if settings.verbose:
                print("BEST FITNESSES:", fitness_list)
                pass
            np.savetxt(settings.output_path + '/fitnesses.dat', fitness_list, delimiter='\n', fmt='%s')


    printIterationInfo(settings, curr_iteration, parent_population, True)
    print('Best GA fitnesses written in fitnesses.dat')

    
def initialise_ga(settings):
    '''

    Initalizes the GA program usign the settings object.
    Picks the corresponding methods for the text strings in the input file that are stored in the settings object from operator_types

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`

    Returns
    -------
    dict
        dictionary containing keys for chosen operations and the method or a list of methods respectively


    '''
    from src.gaapi import operator_types

    chosenSelector, chosenMutator, chosenCrossover, chosenReplacer,chosenGenerator, evaluators = None, None, None, None, None, None
    
    if (settings.selector == operator_types.TOURNAMENT_WOR):
        chosenSelector = selectors.tournamentSelectionWOR

    if (settings.mutator == operator_types.UNI_M):
        chosenMutator = mutators.uniform_mutation
    elif (settings.mutator == operator_types.POLY_M):
        chosenMutator = mutators.polynomial_mutation

    if (settings.crossover == operator_types.UNI_X):
        chosenCrossover = crossovers.uniform_crossover
    elif(settings.crossover == operator_types.SBX):
        chosenCrossover = crossovers.simulated_binary_crossover

    if (settings.replacer == operator_types.ELITIST):
        chosenReplacer = replacers.elitist_replace

    if (settings.replacer== operator_types.NSGA2):
        chosenReplacer = replacers.non_dominated_sorting

    if (settings.generator== operator_types.SWISSPROT):
        chosenGenerator = generators.swissprot_composition_generator
    elif (settings.generator == operator_types.UNBIASED):
        chosenGenerator=generators.unbiased_protein_composition_generator
    else:
        chosenGenerator = generators.UniformCompositionGenerator

    evaluators = []
    for k in settings.evaluators:
        if k in  operator_types.availableEvaluators:
            evaluators.append(operator_types.availableEvaluators[k])
        else:
            print("Could not parse evaluator from setting file")
            sys.exit()

    if len(evaluators) > 1:
        print("More than one evaluator is used. Automatically using non dominated sorting")
        settings.no_evaluators = len(evaluators)
        chosenReplacer = replacers.non_dominated_sorting
        pass
    
    print(settings.composition_residue_indexes)
    print(settings.composition_lower_bounds)
    print(settings.composition_upper_bounds)
    for i,res in enumerate(settings.composition_residue_indexes):
        print("Resindex", res)
        print(i)
        print("First Rotamer",cnts.rotamers[settings.composition_lower_bounds[i]] )
        print("Last Rotamer",cnts.rotamers[settings.composition_upper_bounds[i]] )
        pass

    return {'chosenSelector': chosenSelector, 'chosenMutator': chosenMutator, 'chosenCrossover': chosenCrossover, 'chosenReplacer': chosenReplacer, 'chosenGenerator': chosenGenerator, 'evaluators': evaluators}

    
def printIterationInfo(settings, curr_iteration, pop, ending):
    """

    Prints information about each ga run.

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    curr_iteration: int
        index of current iteration
    pop: list
        contains :class:`src.gaapi.Individual.Individual`
    ending: bool
        to prematurely exit the calculation

    Returns
    -------

    """

    print("CURRENT POPULATION")

    if (settings.backbone_dihedral_optimization or settings.sidechain_dihedral_optimisation):
        if (settings.basilisk_and_sidechains and settings.verbose):
            
            print("ind, ind fitness(es), res, phi, psi, chi angles")
            for i in range (settings.population_size):
                for j in range(0, len(pop[i].dihedral_residue_indexes)):
                    print("{0}, {5}, {1}, {2}, {3}, {4}".format(i, pop[i].dihedral_residue_indexes[j], pop[i].phi_dihedrals[j], pop[i].psi_dihedrals[j], pop[i].chi_angles[j], pop[i].fitnesses))
    
        if (settings.verbose and not settings.basilisk_and_sidechains):
            print("ind, ind fitness(es), res, phi, psi")
            for i in range (settings.population_size):
                for j in range(0, len(pop[i].dihedral_residue_indexes)):
                    print("{0}, {4}, {1}, {2}, {3}".format(i, pop[i].dihedral_residue_indexes[j], pop[i].phi_dihedrals[j], pop[i].psi_dihedrals[j], pop[i].fitnesses))

    print("Current Iteration: ", curr_iteration)
    min_indiv = np.argmin([p.fitnesses[0] for p in pop])
    
    print("Population Average Fitness(es):", np.mean([p.fitnesses for p in pop]))
    print("Population STD(s):", np.std([p.fitnesses for p in pop]))
    if settings.multi_individual:
        print("no best individual for multi-objective runs")
    else:
        print("Best individual:", min_indiv)
    
    if (pop[min_indiv].composition_residue_indexes != None):
        if (settings.multi_individual):
            #print("Best Individual Composition: ", [mi.getResType(pop[min_indiv].mol[0].GetResidue(j)) for j in settings.composition_residue_indexes])
            #print("Best Individual Rotamers: ", [cnts.selected_rotamers[v] for v in pop[min_indiv].composition])
            print("Utopia point individuals")
        else:
            print("Best Individual Composition: ", [mi.getResType(pop[min_indiv].mol.GetResidue(j)) for j in settings.composition_residue_indexes])
            print("Best Individual Rotamers: ", [cnts.selected_rotamers[v] for v in pop[min_indiv].composition])
    
    if (pop[min_indiv].dihedral_residue_indexes != None):
        print('Best ind - PHIs, PSIs')
        print(pop[min_indiv].phi_dihedrals, pop[min_indiv].psi_dihedrals)
        
    print("Best ind - Fitness(es):", pop[min_indiv].fitnesses)

    if settings.write_energies:
        op.writeEnergyLog(path=settings.output_path, energies=pop[min_indiv].energies, iteration=curr_iteration)
            
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdb")

    if (curr_iteration == settings.max_iteration or ending):
        print("OUT: FINAL MOLECULE in {}".format(settings.output_path + '/' + settings.output_file))
        if (settings.multi_individual):
            obConversion.WriteFile(pop[min_indiv].mol[0], settings.output_path + '/' + settings.output_file)
        else:
            obConversion.WriteFile(pop[min_indiv].mol, settings.output_path + '/' + settings.output_file)

    elif (settings.store_intermediate_mol):
        print("OUT: INTERMEDIATE MOLECULE in {}".format(settings.output_path + '/' + "min_mol_iter_" + str(curr_iteration) + ".pdb"))

        if (settings.multi_individual):
            obConversion.WriteFile(pop[min_indiv].mol[0],settings.output_path + '/' + "min_mol_iter_" + str(curr_iteration) + ".pdb")
        else:
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
    #printResidueInfo(settings.initial_molecule)
    print("--- END MOLECULE PRINT --")
    
    print("\nRunning GA")
    mainLoop(settings)
    print("Finished GA in --- {:.2f} seconds ---".format(time.time() - start_time))
