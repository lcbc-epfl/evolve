'''
crossovers.py

@author: Nicholas Browning
'''
import copy
import numpy as np
import Individual
def crossover(settings, individuals, crossover_op, **args):
    
    print("CROSSOVER {}".format(settings.crossover))
    new_population = []
    
    for i in xrange (0, settings.population_size, 2):
        
        if np.random.random() < settings.mating_probability:
            sis, bro = crossover_op(settings, individuals[i], individuals[i + 1], **args)
        else:
            sis, bro = individuals[i], individuals[i + 1]
            
        new_population.append(sis)
        new_population.append(bro)
        
    return new_population
        
def uniform_crossover(settings, mom, dad, **kargs):
        
    bro = Individual.Individual(settings, dad)
    sis = Individual.Individual(settings, mom)
            
            
    bro_phi, bro_psi = bro.phi_dihedrals, bro.psi_dihedrals
    sis_phi, sis_psi = sis.phi_dihedrals, sis.psi_dihedrals
    
        # swap dihedrals
        
    if (settings.dihedral_optimization):
        for i in xrange (0, len(settings.dihedral_residue_indexes)):
                
            if (np.random.random() < settings.genewise_crossover_probability):
                bro.phi_dihedrals[i] = sis_phi[i]
                sis.phi_dihedrals[i] = bro_phi[i]
                
                if (settings.verbose):
                    # TODO
                    pass  
                
            if (np.random.random() < settings.genewise_crossover_probability):
                bro.psi_dihedrals[i] = sis_psi[i]
                sis.psi_dihedrals[i] = bro_psi[i]
                
                if (settings.verbose):
                    # TODO
                    pass
                
    if (settings.composition_optimization):
        
        bro_comp, sis_comp = dad.composition, mom.composition
        
        for i in xrange (0, len(mom.composition)):
            
              if (np.random.random() < settings.genewise_crossover_probability):
                bro.composition[i] = sis_comp[i]
                sis.composition[i] = bro_comp[i]
                
                if (settings.verbose):
                    # TODO
                    pass
            
    return sis, bro


def sim_b(m, d, di, lb, ub):
    if m > d:
        m, d = d, m
        
    if (m == d):
        return m, m
    beta = 1.0 + 2 * min(m - lb, ub - d) / float(d - m)
    alpha = 2.0 - 1.0 / beta ** (di + 1.0)
    u = np.random.random() 
    if u <= (1.0 / alpha):
        beta_q = (u * alpha) ** (1.0 / float(di + 1.0))
    else:
        beta_q = (1.0 / (2.0 - u * alpha)) ** (1.0 / float(di + 1.0))
    bro_val = 0.5 * ((m + d) - beta_q * (d - m))
    bro_val = max(min(bro_val, ub), lb)        
    sis_val = 0.5 * ((m + d) + beta_q * (d - m))
    sis_val = max(min(sis_val, ub), lb)
    if  np.random.random() > 0.5:
        bro_val, sis_val = sis_val, bro_val
    return sis_val, bro_val

def simulated_binary_crossover(settings, mom, dad, **kargs):
    from src import constants
    di = settings.sbx_distribution_index

    bro = Individual.Individual(settings, dad)
    sis = Individual.Individual(settings, mom)
    
    if (settings.dihedral_optimization):
        lower_bound = -180.0
        upper_bound = 180.0
        for i in xrange (0, len(mom.phi_dihedrals)):
        
            m_phi, m_psi = mom.phi_dihedrals[i], mom.psi_dihedrals[i]
            d_phi, d_psi = dad.phi_dihedrals[i], dad.psi_dihedrals[i]
            
            if (np.random.random() < settings.genewise_crossover_probability and m_phi != 999):
            
                sis_phi, bro_phi = sim_b(m_phi, d_phi, di, lower_bound, upper_bound)
                sis.phi_dihedrals[i] = sis_phi
                bro.phi_dihedrals[i] = bro_phi
            
                if (settings.verbose):
                    # TODO
                    pass
             
            if (np.random.random() < settings.genewise_crossover_probability and m_psi != 999):
            
                sis_psi, bro_psi = sim_b(m_psi, d_psi, di, lower_bound, upper_bound)
                sis.psi_dihedrals[i] = sis_psi
                bro.psi_dihedrals[i] = bro_psi
                
                if (settings.verbose):
                    # TODO
                    pass
                
    if (settings.composition_optimization):
        
        lower_comp_bound = settings.composition_lower_bounds
        upper_comp_bound = settings.composition_upper_bounds
        
        print lower_comp_bound, upper_comp_bound
        print mom.composition, dad.composition
        print [constants.rotamers[v] for v in mom.composition], [constants.rotamers[v] for v in dad.composition]
        
        for i in xrange (0, len(mom.composition)):
            
            if (np.random.random() < settings.genewise_crossover_probability):
            
                sis_comp, bro_comp = sim_b(mom.composition[i], dad.composition[i], di, lower_comp_bound[i], upper_comp_bound[i])
                sis.composition[i] = sis_comp
                bro.composition[i] = bro_comp
                
                if (settings.verbose):
                    # TODO
                    pass
            
    
    return sis, bro

    
