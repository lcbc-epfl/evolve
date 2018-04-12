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

def crossover(settings, population, mating_pools):
    
    def sim_b_cross(settings, mom_ang, dad_ang):
        lb = -180
        ub = 180
        di = settings.sbx_distribution_index
        
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
        new_population.append(Individual.Individual(settings, population[i]))


    for i in range (0, len(mating_pools)):
        # loop over residues
        mating_pool = mating_pools[i]
        
        for j in range(settings.population_size, 2):
            
            mom_phi = population[mating_pool[j]].phi_dihedrals[i]
            dad_phi = population[mating_pool[j + 1]].phi_dihedrals[i]
            
            bro_phi, sis_phi = cross(settings, mom_phi, dad_phi)
            
            mom_psi = population[mating_pool[j]].psi_dihedrals[i]
            dad_psi = population[mating_pool[j + 1]].psi_dihedrals[i]
            
            bro_psi, sis_psi = cross(settings, mom_psi, dad_psi)
            
            new_population[j].phi_dihedrals[i] = bro_phi
            new_population[j + 1].phi_dihedrals[i] = sis_phi
            
            new_population[j].psi_dihedrals[i] = bro_psi
            new_population[j + 1].psi_dihedrals[i] = sis_psi
            
    
    return new_population

def evaluate(settings, individuals):
    evals.amber_energy_simplified(settings, individuals, 0, pop_start=0)
   
    for i in xrange(0, len(individuals)):
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat("pdb")
        obconv.WriteFile(individuals[i].mol, directory + "/min_struct.pdb")
        op.runtleap(work_dir=directory + "/", mol_file='min_struct.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in", inpcrd_out="min_struct.inpcrd", prmtop_out="min_struct.prmtop")
        op.runMMPBSA(work_dir=directory + "/", prmtop="min_struct.prmtop", inpcrd="min_struct.inpcrd")
    
        decomp_dict = op.parseMMBPSA_total_decomp("mmpbsa_decomp.dat", individuals[i].mol.NumResidues())
        
        for j in range(0, len(settings.dihedral_residue_indexes)):
            individuals[i].decomp_fitnesses[j] = return_dict[j][5]
        

def replace (settings, parents, children):
    return children

def mutate (settings, population):
    return population

def tournament_select (settings, population, residue_index):
    
    fitnesses = [population[i].decomp_fitnesses[residue_index] for i in range(settings.population_size)]
    i = 0
    mating_pool = np.zeros(settings.population_size)
    
    while i < settings.population_size:
        candidates = np.random.choice(settings.population_size, settings.tournament_size)
        
        best_candidate = -1
        
        for i in range(settings.tournament_size):
            if (best_candidate == -1 or fitnesses[candidates[i]] < fitnesses[best_candidate]):
                best_candidate = candidates[i]
        
        mating_pool[i] = best_candidate
    
        i += 1
        
    
    return mating_pool


def printPop(settings, population):
    for i in range (0, len(population)):
        print i, population[i].fitnesses, "DECOMP:" , population[i].decomp_fitnesses
        
    best_idx = np.argmin([indiv.fitnesses[0] for indiv in population])
    print "BEST:", population[best_idx].fitnesses[0], "DECOMP:", population[best_idx].decomp_fitnesses
        
def mainLoop(settings):
    
    parents = []
    
    for i in xrange (settings.population_size):
        
       indiv = Individual.Individual(settings)
       
       indiv.decomp_fitnesses = np.zeros(len(settings.dihedral_residue_indexes))
       
       parents.append(indiv)
       
    generators.BasiliskSideChainDihedralsGenerator(settings, parents)
    
    evaluate(settings, parents)
    
    curr_iteration = 0
    
    while (curr_iteration < settings.max_iteration):
        
        # loop over residues
        mating_pools = []
        for i in range(len(settings.dihedral_residue_indexes)):
            mating_pools.append(tournament_select(settings, parents, i))
            
        children = crossover(settings, parents, mating_pools)
        children = mutate(settings, children)
            
        evaluate(settings, children)
        printPop(settings, children)
        
        parent_population = replace(settings, parents, children)
        
        
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
