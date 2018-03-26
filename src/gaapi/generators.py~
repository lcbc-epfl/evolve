'''
crossovers.py

@author: Nicholas Browning
'''
import Individual
import constants
import numpy as np

def initialisePopulation(settings):
    
    individuals = []
    
    for i in xrange (0, settings.population_size):
        indiv = Individual.Individual(settings)
        
        individuals.append(indiv)
        
    return individuals
        
def defaultGenerator(settings):
    pass

def MonteCarloDihedralGenerator(settings, individuals):
    
#     [MC_GENERATE]
# mc_generate_dihedrals = True
# distribution_path = 's'
# mcmove_lbound = -45.0
# mcmove_ubound = 45.0
# num_mc_steps = 1000
# dihedral_probability_pointers = 1 2 3 1

    from src.montecarlo import dihedral_opt as mcopt
 
    probs = mcopt.loadDefaultProbabilityDistributions()
    
    num_phipsi = len(settings.dihedral_residue_indexes)
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi * settings.population_size, 2))
    
    print "INITIAL DIHEDRALS"
    
    for i in xrange (settings.population_size):
        print initial_solution[i * num_phipsi: (i + 1) * num_phipsi]
    move_lows = [settings.mcmove_lbound] * num_phipsi * settings.population_size
    move_highs = [settings.mcmove_ubound] * num_phipsi * settings.population_size
    
    # TODO do this once with N*M array or per individual? 
    # Can pack intial solution as an N*M array without problems...probably easiest - need to also create temp prob_pointer array of same size
    
    pointers = []
    
    for i in xrange (0, settings.population_size):
        for j in xrange (0, len (settings.dihedral_probability_pointers)):
            pointers.append(settings.dihedral_probability_pointers[j])
    
    dihedrals = mcopt.optimisePsiPhi(probs, pointers, initial_solution, settings.num_mc_steps, move_lows, move_highs)
    print "MC_OPT DIHEDRALS"
    
    for i in xrange (settings.population_size):
        print dihedrals[i * num_phipsi: (i + 1) * num_phipsi]
        
    for i in xrange (0, settings.population_size):
        indiv = individuals[i]
        for j in xrange (0, len(indiv.phi_dihedrals)):
            indiv.phi_dihedrals[j] = dihedrals[i * num_phipsi + j][0]
            indiv.psi_dihedrals[j] = dihedrals[i * num_phipsi + j][1]
        
def UniformDihedralGenerator(settings, individuals):
    
    num_phipsi = len(settings.dihedral_residue_indexes)
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi, 2))
    
    for i in xrange (0, settings.population_size):
        indiv = individuals[i]
        for j in xrange (0, len(indiv.phi_psi_dihedrals)):
            indiv.phi_psi_dihedrals[j] = (dihedrals[j][0], dihedrals[j][1])
        indiv.applyPhiPsiDihedrals()
    

def UniformCompositionGenerator(settings, individuals):
    
    for i in xrange (0, settings.population_size):
        indiv = individuals[i]
        for j in xrange (0, len(settings.composition_residue_indexes)):
            indiv.composition[j] = np.random.randint(low=settings.composition_lower_bounds[j], high=settings.composition_upper_bounds[j])
        
        
    
            
        
        
