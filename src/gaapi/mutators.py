'''
crossovers.py

@author: Nicholas Browning
'''

import numpy as np
import copy
import Individual

def mutate(settings, individuals, mutator_op, **args):
    mutated_pop = []
    
    for i in xrange (0, settings.population_size):
        if (np.random.random() < settings.mutation_probability):
            mutated_pop.append(mutator_op(settings, individuals[i], **args))
        else:
            mutated_pop.append(individuals[i])
            
    return mutated_pop
    

def poly(indiv_i, xl, xu, eta):
    x = indiv_i
    delta_1 = (x - xl) / (xu - xl)
    delta_2 = (xu - x) / (xu - xl)
    rand = np.random.random()
    
    mut_pow = 1.0 / (eta + 1.)
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
    
def polynomial_dihedral(settings, individual, **args):
 
    eta = settings.poly_eta
    mut_rate = settings.genewise_mutation_probability

    mutant = Individual.Individual(settings, individual)
    
    if (settings.dihedral_optimization):
        
        lower_bound = -180
        upper_bound = 180
    
        for i in xrange (0, len(settings.dihedral_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                psi_m = poly(individual.psi_dihedrals[i], lower_bound, upper_bound, eta)
                phi_m = poly(individual.phi_dihedrals[i], lower_bound, upper_bound, eta)

                if (settings.verbose):
                    print "DHDRL: Mutating!", individual.psi_dihedrals[i], psi_m, individual.phi_dihedrals[i], phi_m
                individual.psi_dihedrals[i] = psi_m

                individual.phi_dihedrals[i] = phi_m
                
    if (settings.composition_optimization):
        
        lower_comp_bound = settings.composition_lower_bounds
        upper_comp_bound = settings.composition_upper_bounds
    
        for i in xrange (0, len(settings.composition_residue_indexes)):
        
            if np.random.random() <= mut_rate:
            
                comp_m = poly(individual.composition[i], lower_comp_bound[i], upper_comp_bound[i], eta)

                if (settings.verbose):
                    print "CMPSTN: Mutating!", individual.composition[i], comp_m
                individual.composition[i] = comp_m

    return mutant
