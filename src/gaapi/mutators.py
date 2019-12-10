'''
crossovers.py

@author: Nicholas Browning
'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import range
from past.utils import old_div
import numpy as np
import copy
from . import Individual

def mutate(settings, individuals, mutator_op, **args):
    print("MUTATOR {}".format(settings.mutator))
    mutated_pop = []
    
    for i in range (0, settings.population_size):
        if (np.random.random() < settings.mutation_probability):
            if (settings.verbose):
                print("MUTATE Individual: ", i)
            mutated_pop.append(mutator_op(settings, individuals[i], **args))
        else:
            mutated_pop.append(individuals[i])
            
    return mutated_pop

def uniform_mutation(settings, individual):

    mut_rate = settings.genewise_mutation_probability

    mutant = Individual.Individual(settings, individual)
    
    if (settings.backbone_dihedral_optimization):
        
        lower_bound = -180
        upper_bound = 180
    
        for i in range (0, len(settings.dihedral_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                psi_m = np.random.rand() * 360 + lower_bound
                phi_m = np.random.rand() * 360 + lower_bound

                if (settings.verbose):
                    print("DHDRL: Mutating!", mutant.psi_dihedrals[i], psi_m, mutant.phi_dihedrals[i], phi_m)
                mutant.psi_dihedrals[i] = psi_m
                mutant.phi_dihedrals[i] = phi_m
                
    if (settings.composition_optimization):
        
        lower_comp_bound = settings.composition_lower_bounds
        upper_comp_bound = settings.composition_upper_bounds
    
        for i in range (0, len(settings.composition_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                comp_m = np.random.randint(low=lower_comp_bound[i], high=upper_comp_bound[i])

                if (settings.verbose):
                    print("Position (", i, " ): " , mutant.composition[i], "->", comp_m)
                    
                mutant.composition[i] = comp_m

    return mutant
    

def poly(indiv_i, xl, xu, eta):
    
    x = float(indiv_i)
    xu = float(xu)
    xl = float(xl)
    eta = float(eta)
    
    delta_1 = old_div((x - xl), (xu - xl))
    delta_2 = old_div((xu - x), (xu - xl))
    rand = np.random.random()
    
    mut_pow = old_div(1.0, (eta + 1.))
    if rand < 0.5:
        xy = 1.0 - delta_1
        val = 2.0 * rand + (1.0 - 2.0 * rand) * xy ** (eta + 1)
        delta_q = val ** mut_pow - 1.0
    else:
        xy = 1.0 - delta_2
        val = 2.0 * (1.0 - rand) + 2.0 * (rand - 0.5) * xy ** (eta + 1)
        delta_q = 1.0 - val ** mut_pow

    x = x + delta_q * (xu - xl)
    x = min(max(x, xl), xu)
    return x
    
def polynomial_mutation(settings, individual, **args):
 
    eta = settings.poly_eta
    mut_rate = settings.genewise_mutation_probability

    mutant = Individual.Individual(settings, individual)
    
    if (settings.backbone_dihedral_optimization):
        
        lower_bound = -180
        upper_bound = 180
    
        for i in range (0, len(settings.dihedral_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                psi_m = poly(individual.psi_dihedrals[i], lower_bound, upper_bound, eta)
                phi_m = poly(individual.phi_dihedrals[i], lower_bound, upper_bound, eta)

                if (settings.verbose):
                    print("DHDRL: Mutating!", individual.psi_dihedrals[i], psi_m, individual.phi_dihedrals[i], phi_m)
                mutant.psi_dihedrals[i] = psi_m

                mutant.phi_dihedrals[i] = phi_m
                
    if (settings.composition_optimization):
        
        lower_comp_bound = settings.composition_lower_bounds
        upper_comp_bound = settings.composition_upper_bounds
    
        for i in range (0, len(settings.composition_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                comp_m = poly(individual.composition[i], lower_comp_bound[i], upper_comp_bound[i], eta)
                
                if (settings.verbose):
                    print("CMPSTN: Mutating!", individual.composition[i], comp_m)
                    
                mutant.composition[i] = int(comp_m)

    return mutant
