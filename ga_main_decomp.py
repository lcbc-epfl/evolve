'''
main.py

@author: Nicholas Browning
'''

import numpy as np
import argparse
import openbabel
import time
import json

import src.outprocesses as op

from src import JobConfig
from src.gaapi import Individual
from src.gaapi import selectors
from src.gaapi import mutators
from src.gaapi import crossovers
from src.gaapi import generators
from src.gaapi import replacers
import evaluators as evals
import printResidueInfo as pri

def replace(settings, parents, children, **args):

    elitist_factor = settings.elitist_factor
    num_elites = int(elitist_factor * settings.population_size)

    print("ELITIST REPLACER with {} elite parents".format(num_elites))
    
    parent_fitnesses = [indiv.fitnesses[0] for indiv in parents]
    child_fitnesses = [indiv.fitnesses[0] for indiv in children]
    
    sorted_parent_indexes = np.argsort(parent_fitnesses)
    sorted_child_indexes = np.argsort(child_fitnesses)
    
    if (settings.verbose):
        print "elite parent indexes:", sorted_parent_indexes[:num_elites]
        print "elite children indexes:", sorted_child_indexes[:settings.population_size - num_elites]
    
    elitist_population = []
    
    for i in xrange (num_elites):
        indiv = Individual.Individual(settings, parents[sorted_parent_indexes[i]])
        indiv.decomp_fitnesses = parents[sorted_parent_indexes[i]].decomp_fitnesses
        elitist_population.append(indiv)
    
    for i in xrange (settings.population_size - num_elites):
        indiv = Individual.Individual(settings, children[sorted_child_indexes[i]])
        indiv.decomp_fitnesses = children[sorted_child_indexes[i]].decomp_fitnesses
        elitist_population.append(indiv)
        
        
    return elitist_population
        
def uniform_point_mutation(settings, individual):

    mut_rate = settings.genewise_mutation_probability

    mutant = Individual.Individual(settings, individual)
    mutant.decomp_fitnesses = individual.decomp_fitnesses
     
    lower_bound = -180
    upper_bound = 180
    
    for i in xrange (0, len(settings.dihedral_residue_indexes)):
        
        if np.random.random() <= mut_rate:
            
            psi_m = np.random.rand() * 360 + lower_bound
            phi_m = np.random.rand() * 360 + lower_bound

            mutant.psi_dihedrals[i] = psi_m
            mutant.phi_dihedrals[i] = phi_m                

    return mutant

def helix_unit_mutation(settings, individual):
    from math import degrees
    
    mutant = Individual.Individual(settings, individual)
    mutant.decomp_fitnesses = individual.decomp_fitnesses
     
    chis, bbs, lls = generators.getBasiliskSample(mutant.mol)
    
    middle = np.random.randint(low=1, high=len(settings.dihedral_residue_indexes) - 1)
    print "3-point helix mutation at: ", middle
    for i in range (middle - 1, middle + 2):
        mutant.phi_dihedrals[i] = degrees(bbs[i][0]) - 180
        mutant.psi_dihedrals[i] = degrees(bbs[i][1]) - 180
        
    return mutant
   
        
def crossover(settings, population, mating_pools):
    
    def sim_b_cross(settings, mom_ang, dad_ang):
        lb = -180
        ub = 180 
        di = 4.0
        
        if (np.random.rand() > settings.genewise_crossover_probability):
            return mom_ang, dad_ang
        
        if mom_ang > dad_ang:
            mom_ang, dad_ang = dad_ang, mom_ang
            
        if (mom_ang == dad_ang):
            return mom_ang, dad_ang
        
        beta = 1.0 + 2 * min(mom_ang - lb, ub - dad_ang) / float(dad_ang - mom_ang)
        alpha = 2.0 - 1.0 / beta ** (di + 1.0)
        u = np.random.random() 
        if u <= (1.0 / alpha):
            beta_q = (u * alpha) ** (1.0 / float(di + 1.0))
        else:
            beta_q = (1.0 / (2.0 - u * alpha)) ** (1.0 / float(di + 1.0))
        bro_val = 0.5 * ((mom_ang + dad_ang) - beta_q * (dad_ang - mom_ang))
        bro_val = max(min(bro_val, ub), lb)        
        sis_val = 0.5 * ((mom_ang + dad_ang) + beta_q * (dad_ang - mom_ang))
        sis_val = max(min(sis_val, ub), lb)
        if  np.random.random() > 0.5:
            bro_val, sis_val = sis_val, bro_val
        return sis_val, bro_val


    def uni_cross(settings, mom_ang, dad_ang):
        if (np.random.rand() < settings.genewise_crossover_probability):
            return dad_ang, mom_ang
        else:
            return mom_ang, dad_ang
        
    cross = sim_b_cross
            
    new_population = []
    
    for i in range (0, settings.population_size):
        new_indiv = Individual.Individual(settings, population[i])
        
        new_indiv.decomp_fitnesses = population[i].decomp_fitnesses
        
        new_population.append(new_indiv)


    for i in range (0, len(mating_pools)):
        # loop over residues
        mating_pool = mating_pools[i]
        # print "mat pool:", i , mating_pool
        
        for j in range(0, settings.population_size, 2):
            
            mom_phi = population[mating_pool[j]].phi_dihedrals[i]
            dad_phi = population[mating_pool[j + 1]].phi_dihedrals[i]
            
            bro_phi, sis_phi = cross(settings, mom_phi, dad_phi)
            
            mom_psi = population[mating_pool[j]].psi_dihedrals[i]
            dad_psi = population[mating_pool[j + 1]].psi_dihedrals[i]
            
            bro_psi, sis_psi = cross(settings, mom_psi, dad_psi)
            
            # print mom_phi, "->", sis_phi 
            
            new_population[j].phi_dihedrals[i] = bro_phi
            new_population[j + 1].phi_dihedrals[i] = sis_phi
            
            new_population[j].psi_dihedrals[i] = bro_psi
            new_population[j + 1].psi_dihedrals[i] = sis_psi
            
    
    return new_population

def evaluate(settings, individuals):
    evals.amber_energy_simplified(settings, individuals, 0, pop_start=0)
    directory = settings.output_path + "/amber_run"
    
    for i in xrange(0, len(individuals)):
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat("pdb")
        obconv.WriteFile(individuals[i].mol, directory + "/min_struct.pdb")
        op.runtleap(work_dir=directory + "/", mol_file='min_struct.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in", inpcrd_out="min_struct.inpcrd", prmtop_out="min_struct.prmtop")
        op.runMMPBSA(work_dir=directory + "/", prmtop="min_struct.prmtop", inpcrd="min_struct.inpcrd")
    
        decomp_dict = op.parseMMBPSA_total_decomp(directory + "/mmpbsa_decomp.dat", individuals[i].mol.NumResidues())
        
        for j in range(0, len(settings.dihedral_residue_indexes)):
            if (decomp_dict is None):
                individuals[i].decomp_fitnesses[j] = 9999
            else:
                individuals[i].decomp_fitnesses[j] = decomp_dict[j][5]


def mutate (settings, population):
    mutated_pop = []
    
    for i in xrange (0, settings.population_size):
        if (np.random.random() < settings.mutation_probability):
            mutated_pop.append(uniform_point_mutation(settings, population[i]))
#             if (np.random.random() < 0.5):
#                 mutated_pop.append(uniform_point_mutation(settings, population[i]))
#             else:
#                 mutated_pop.append(helix_unit_mutation(settings, population[i]))
        else:
            mutated_pop.append(population[i])
            
    return mutated_pop

def tournament_select (settings, population, residue_index):
    
    fitnesses = [population[i].decomp_fitnesses[residue_index] for i in range(settings.population_size)]
    
    mating_pool = np.zeros(settings.population_size, dtype=np.int32)
    i = 0
    while i < settings.population_size:
 
        candidates = np.random.choice(settings.population_size, settings.tournament_size)
        
        best_candidate = -1
        
        for j in range(settings.tournament_size):
            if (best_candidate == -1 or fitnesses[candidates[j]] < fitnesses[best_candidate]):
                best_candidate = candidates[j]
        
        mating_pool[i] = best_candidate
    
        i += 1
        
    
    return mating_pool


def printPop(settings, population, curr_iter):
    for i in range (0, len(population)):
        print i, population[i].fitnesses, "DECOMP:" , population[i].decomp_fitnesses
        
    best_idx = np.argmin([indiv.fitnesses[0] for indiv in population])
    print "Iteration", curr_iter
    print "BEST:", population[best_idx].fitnesses[0], "DECOMP:", population[best_idx].decomp_fitnesses
    
    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    
    obconv.WriteFile(population[best_idx].mol, settings.output_path + '/' + "min_mol_iter_" + str(curr_iter) + ".pdb")

        
def mainLoop(settings):
    
    initial_indiv = Individual.Individual(settings)
    initial_indiv.decomp_fitnesses = np.zeros(len(settings.dihedral_residue_indexes))
   
    evaluate(settings, [initial_indiv])
    with open(settings.output_path + '/initial_fitness.dat', 'w') as outfile:  
        json.dump(initial_indiv.fitnesses[0], outfile)
        
    parents = []
    
    for i in xrange (settings.population_size):
        
       indiv = Individual.Individual(settings)
       
       indiv.decomp_fitnesses = np.zeros(len(settings.dihedral_residue_indexes))
       
       parents.append(indiv)
    
    # Generate dihedrals by Monte Carlo
    generators.MonteCarloDihedralGenerator(settings, parents)
    # generators.BasiliskSideChainDihedralsGenerator(settings, parents)
    
    evaluate(settings, parents)
    
    print "INITIAL_POPULATION"
    printPop(settings, parents, -1)
    
    curr_iteration = 0
    
    while (curr_iteration < settings.max_iteration):
        
        # loop over residues
        mating_pools = []
        for i in range(len(settings.dihedral_residue_indexes)):
            mating_pools.append(tournament_select(settings, parents, i))
            
        children = crossover(settings, parents, mating_pools)
        children = mutate(settings, children)
            
        evaluate(settings, children)
        
        parents = replace(settings, parents, children)
        
        printPop(settings, parents, curr_iteration)
        
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
    
    print("\nRunning GA")
    mainLoop(settings)
    print("Finished GA in --- {:.2f} seconds ---".format(time.time() - start_time))
