'''

not sure what this dois?

.. codeauthor:: Nicholas Browning
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

import copy

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
from src import MoleculeInfo as mi

def replace(settings, parents, children, **args):

    elitist_factor = settings.elitist_factor
    num_elites = int(elitist_factor * settings.population_size)

    print("ELITIST REPLACER with {} elite parents".format(num_elites))
    
    parent_fitnesses = np.asarray([indiv.fitnesses[0] for indiv in parents])
    child_fitnesses = np.asarray([indiv.fitnesses[0] for indiv in children])
    
    sorted_parent_indexes = np.argsort(parent_fitnesses)
    
    print(parent_fitnesses[sorted_parent_indexes])
    
    sorted_child_indexes = np.argsort(child_fitnesses)
    
    print(child_fitnesses[sorted_child_indexes])
    
    if (settings.verbose):
        print("elite parent indexes:", sorted_parent_indexes[:num_elites])
        print("elite children indexes:", sorted_child_indexes[:settings.population_size - num_elites])
    
    elitist_population = []
    
    for i in range (num_elites):
       # indiv = Individual.Individual(settings, parents[sorted_parent_indexes[i]])
        # indiv.total_decomp_fitness = parents[sorted_parent_indexes[i]].total_decomp_fitness
        elitist_population.append(parents[sorted_parent_indexes[i]])
    
    for i in range (settings.population_size - num_elites):
        # indiv = Individual.Individual(settings, children[sorted_child_indexes[i]])
        # indiv.total_decomp_fitness = children[sorted_child_indexes[i]].total_decomp_fitness
        elitist_population.append(children[sorted_child_indexes[i]])
        
        
    return elitist_population

def mutate (settings, population):
    mutated_pop = []
    
    for i in range (0, settings.population_size):
        if (np.random.random() < settings.mutation_probability):
            mutated_pop.append(uniform_point_mutation(settings, population[i]))
#             if (np.random.random() < 0.5):
#                 mutated_pop.append(uniform_point_mutation(settings, population[i]))
#             else:
#                 mutated_pop.append(helix_unit_mutation(settings, population[i]))
        else:
            mutated_pop.append(population[i])
            
    return mutated_pop
       
def uniform_point_mutation(settings, individual):

    mut_rate = settings.genewise_mutation_probability

    mutant = Individual.Individual(settings, individual)
    mutant.total_decomp_fitness = copy.deepcopy(individual.total_decomp_fitness)
     
    lower_bound = -180
    upper_bound = 180
    
    for i in range (0, len(settings.dihedral_residue_indexes)):
        
        if np.random.random() <= mut_rate:
            
            psi_m = lower_bound + np.random.rand() * 360
            phi_m = lower_bound + np.random.rand() * 360 

            if (mutant.psi_dihedrals[i] is not None): 
                mutant.psi_dihedrals[i] = psi_m
                
            if (mutant.phi_dihedrals[i] is not None): 
                mutant.phi_dihedrals[i] = phi_m
        
        if np.random.random() < mut_rate:
            for j in range (len(mutant.chi_angles[i])):
                mut_rate_chi = old_div(0.5, len(mutant.chi_angles[i]))
                if np.random.random() < mut_rate_chi:
                    mutant.chi_angles[i][j] = lower_bound + np.random.rand() * 360

    return mutant

def helix_unit_mutation(settings, individual):
    from math import degrees
    
    mutant = Individual.Individual(settings, individual)
    mutant.total_decomp_fitness = copy.deepcopy(individual.total_decomp_fitness)
     
    chis, bbs, lls = generators.getBasiliskSample(mutant.mol)
    
    middle = np.random.randint(low=1, high=len(settings.dihedral_residue_indexes) - 1)
    print("3-point helix mutation at: ", middle)
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
        
        beta = 1.0 + 2 * min(mom_ang - lb, ub - dad_ang) / (dad_ang - mom_ang)
        alpha = 2.0 - old_div(1.0, beta ** (di + 1.0))
        u = np.random.random() 
        if u <= (old_div(1.0, alpha)):
            beta_q = (u * alpha) ** (old_div(1.0, float(di + 1.0)))
        else:
            beta_q = (old_div(1.0, (2.0 - u * alpha))) ** (old_div(1.0, (di + 1.0)))
        bro_val = 0.5 * ((mom_ang + dad_ang) - beta_q * (dad_ang - mom_ang))
        bro_val = max(min(bro_val, ub), lb)        
        sis_val = 0.5 * ((mom_ang + dad_ang) + beta_q * (dad_ang - mom_ang))
        sis_val = max(min(sis_val, ub), lb)
        
        # if  np.random.random() > 0.5:
        #    bro_val, sis_val = sis_val, bro_val
            
        if (sis_val < -180):
            sis_val = -180.0
        if (sis_val > 180):
            sis_val = 180.0
        
        if (bro_val < -180):
            bro_val = -180.0
        if (bro_val > 180):
            bro_val = 180.0
            
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
        
        new_indiv.total_decomp_fitness = copy.deepcopy(population[i].total_decomp_fitness)
        
        new_population.append(new_indiv)


    for i in range (0, len(mating_pools)):
        # loop over residues
        mating_pool = mating_pools[i]
        # print "mat pool:", i , mating_pool
        
        for j in range(0, settings.population_size, 2):
            
            # cross backbone phi-psi
            
            mom_phi = population[mating_pool[j]].phi_dihedrals[i]
            dad_phi = population[mating_pool[j + 1]].phi_dihedrals[i]
            
            if (mom_phi is not None and  dad_phi is not None):
                bro_phi, sis_phi = cross(settings, mom_phi, dad_phi)
                new_population[j].phi_dihedrals[i] = bro_phi
                new_population[j + 1].phi_dihedrals[i] = sis_phi
            
            mom_psi = population[mating_pool[j]].psi_dihedrals[i]
            dad_psi = population[mating_pool[j + 1]].psi_dihedrals[i]
            
            if (mom_psi is not None and  dad_psi is not None):
                bro_psi, sis_psi = cross(settings, mom_psi, dad_psi)
                new_population[j].psi_dihedrals[i] = bro_psi
                new_population[j + 1].psi_dihedrals[i] = sis_psi
            
            # cross sidechain
            
            mom_chis = population[mating_pool[j]].chi_angles[i]
            dad_chis = population[mating_pool[j + 1]].chi_angles[i]
            
            if len(mom_chis) == 0:
                continue
            
            for k in range (len(mom_chis)):
                bro_chi, sis_chi = cross(settings, mom_chis[k], dad_chis[k])
                new_population[j].chi_angles[i][k] = sis_chi
                new_population[j + 1].chi_angles[i][k] = bro_chi
            
            
    
    return new_population

def evaluate(settings, individuals):
    import os
    
    evals.amber_energy_simplified(settings, individuals, 0, pop_start=0)
    
    directory = settings.output_path + "/amber_run"
    
    for i in range(0, len(individuals)):
        print("Computing MMPBSA", i)
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat("pdb")
        obconv.WriteFile(individuals[i].mol, directory + "/min_struct.pdb")
        
        op.runtleap(work_dir=directory + "/", mol_file='min_struct.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in", inpcrd_out="min_struct.inpcrd", prmtop_out="min_struct.prmtop")
        
        op.runMMPBSA(work_dir=directory + "/", prmtop="min_struct.prmtop", inpcrd="min_struct.inpcrd")
        
        pair_decomp, decomp_dict = op.parseMMPBSA_pairwise_total(directory + "/mmpbsa_decomp.dat", individuals[i].mol.NumResidues())
        
        total_e = op.parseMMPBSA(directory + "/mmpbsa.dat")
        
        for j in range(0, len(settings.dihedral_residue_indexes)):
            if (decomp_dict is None):
                individuals[i].total_decomp_fitness[j] = 9999
            else:
                individuals[i].total_decomp_fitness[j] = np.sum(pair_decomp[j]) - pair_decomp[j, j]
        
        individuals[i].setFitness(0, total_e)
        
        print(individuals[i].fitnesses)
        print(individuals[i].total_decomp_fitness)
        
        if os.path.isfile(directory + "/mmpbsa_decomp.dat"):
            os.remove(directory + "/mmpbsa_decomp.dat")
            
        if os.path.isfile(directory + "/mmpbsa.dat"):
           os.remove(directory + "/mmpbsa.dat")
            
        if os.path.isfile(directory + "/mmpbsa.log"):
           os.remove(directory + "/mmpbsa.log")


def tournament_select (settings, population, residue_index):
    
    fitnesses = [population[i].total_decomp_fitness[residue_index] for i in range(settings.population_size)]
    
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
        print(i, population[i].fitnesses, "DECOMP:" , population[i].total_decomp_fitness)
        
    best_idx = np.argmin([indiv.fitnesses[0] for indiv in population])
    print("Iteration", curr_iter)
    print("BEST:", population[best_idx].fitnesses[0], "DECOMP:", population[best_idx].total_decomp_fitness)
    
    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    
    obconv.WriteFile(population[best_idx].mol, settings.output_path + '/' + "min_mol_iter_" + str(curr_iter) + ".pdb")
        
def mainLoop(settings):
    
    settings.chi_dihedral_atom_idxs = mi.getChiDihedralAtomIndexes(settings.initial_molecule, settings.dihedral_residue_indexes)
    settings.backbone_dihedral_atom_idxs = mi.getPhiPsiDihedralAtomIndexes(settings.initial_molecule, settings.dihedral_residue_indexes)
    
    print("Backbone dihedral atom idxs:")
    print(settings.backbone_dihedral_atom_idxs)
    print("Chi dihedral atom idxs:")
    print(settings.chi_dihedral_atom_idxs)
    
    initial_indiv = Individual.Individual(settings)
    initial_indiv.total_decomp_fitness = np.zeros(len(settings.dihedral_residue_indexes))
    evaluate(settings, [initial_indiv])
    
    with open(settings.output_path + '/initial_fitness.dat', 'w') as outfile:  
        json.dump(initial_indiv.fitnesses[0], outfile)
        
    parents = []
    
    for i in range (settings.population_size):
        
       indiv = Individual.Individual(settings)
       
       indiv.total_decomp_fitness = np.zeros(len(settings.dihedral_residue_indexes))
       
       parents.append(indiv)
    
    # Generate dihedrals by Monte Carlo
    print("Generating backbone via MC")
    generators.MonteCarloDihedralGenerator(settings, parents)
    print("Generating sidechain via Basilisk")
    generators.BasiliskSideChainDihedralsGenerator(settings, parents)
    
    for i in range (settings.population_size):
        print("setting initial population", i)
        parents[i].applyChiDihedrals(settings)
        parents[i].applyPhiPsiDihedrals(settings)
            
    evaluate(settings, parents)
    
    print("INITIAL_POPULATION")
    printPop(settings, parents, -1)
    
    curr_iteration = 0
    
    while (curr_iteration < settings.max_iteration):
        
        # loop over residues
        mating_pools = []
        for i in range(len(settings.dihedral_residue_indexes)):
            mating_pools.append(tournament_select(settings, parents, i))
            
        children = crossover(settings, parents, mating_pools)
        children = mutate(settings, children)
        
        for i in range (settings.population_size):
            children[i].applyPhiPsiDihedrals(settings)
            children[i].applyChiDihedrals(settings)
        
        evaluate(settings, children)
    
        print("Updating angles, dihedrals")
        for i in range (settings.population_size):
            children[i].updateChiAngles(settings)
            children[i].updatePhiPsiDihedrals(settings)

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
