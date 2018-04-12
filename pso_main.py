'''
main.py

@author: Nicholas Browning
'''

import numpy as np
import argparse
import openbabel
import time
import json

from src.gaapi import generators
from src import JobConfig
import evaluators as evals
import printResidueInfo as pri
import evaluators as evals

from src.gaapi import Individual


def mainLoop(settings):
    
    parent_population = []
   
    for i in xrange (settings.population_size):
       indiv = Individual.PSOIndividual(settings)
       
       # initialise velolcities
       indiv.phi_velocity_vector = 0.1 * (np.random.rand(len(settings.dihedral_residue_indexes)) * 360) - 180
       indiv.phi_velocity_vector = 0.1 * (np.random.rand(len(settings.dihedral_residue_indexes)) * 360) - 180
       
       # initialise topological neighbours
       
       indexes = np.arange(settings.population_size)
       indexes = np.delete(indexes, i)
       np.random.shuffle(indexes)
       num_neighbours = np.random.randint(low=3, high=7)
       indiv.topological_neighbours = indexes[:num_neighbours]
       indiv.topological_best_neighbour = indiv.topological_neighbours[0]
       
       parent_population.append(indiv)
       
    generators.BasiliskSideChainDihedralsGenerator(settings, parent_population)
   
    curr_iteration = 0
   
    best_fitnesses = [np.min([p.fitnesses[0] for p in parent_population])]
    
    evals.amber_energy_simplified(settings, parent_population, 0)
   
    while (curr_iteration < settings.max_iteration):
        
        print ('Current Iteration:', curr_iteration)
        
        for i in xrange (settings.population_size):
            parent_population[i].updatePersonalBest()
            parent_population[i].updateTopologicalBestNeighbour(parent_population)
            parent_population[i].updateVelocity(parent_population, 0.6, 0.5, 0.5)
            print parent_population[i].phi_velocity_vector
            print parent_population[i].psi_velocity_vector
            parent_population[i].updatePosition()
    
        evals.amber_energy_simplified(settings, parent_population, 0)
       
        curr_iteration += 1
       

    
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
    
    print("\nRunning PSO")
    mainLoop(settings)
    print("Finished PSO in --- {:.2f} seconds ---".format(time.time() - start_time))
