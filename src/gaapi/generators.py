'''
generators.py

@author: Nicholas Browning
'''

# Import modules
import Individual
import constants
import numpy as np

from src import MoleculeInfo as mi
import openbabel

# Initialisation of the first population
def initialisePopulation(settings):
    
    individuals = []
    
    for i in xrange (0, settings.population_size):
        indiv = Individual.Individual(settings)
        
        individuals.append(indiv)
        
    return individuals

# TO DO
def defaultGenerator(settings):
    pass

# Generate dihedrals with uniform distribution
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

            
# Generate dihedrals by Monte Carlo
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
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi * len(individuals), 2))
    
    print "INITIAL DIHEDRALS UNIFORMLY AND RANDOMLY GENERATED"

    if (settings.verbose):
        for i in xrange(0, len(individuals)):
            print initial_solution[i * num_phipsi: (i + 1) * num_phipsi]
            
    move_lows = [settings.mcmove_lbound] * num_phipsi * settings.population_size
    move_highs = [settings.mcmove_ubound] * num_phipsi * settings.population_size
    
    # TODO do this once with N*M array or per individual? 
    # Can pack intial solution as an N*M array without problems...probably easiest - need to also create temp prob_pointer array of same size
    
    pointers = []

    # WOULD NEED TO CHANGE THIS FUNCTION IN ORDER NOT TO SPECIFY THE POINTERS IN THE SETTING FILE, WOULD HAVE AN AUTOMATICE PROCESS COMING FROM guessResidueDIHProbPointers.py
    for i in xrange(0, len(individuals)):
        for j in xrange(0, len(settings.dihedral_probability_pointers)):
            pointers.append(settings.dihedral_probability_pointers[j])
    
    dihedrals = mcopt.optimisePsiPhi(probs, pointers, initial_solution, settings.num_mc_steps, move_lows, move_highs)
    print "MC_OPTIMIZED DIHEDRALS GENERATED"

    if (settings.verbose):
        for i in xrange (0, len(individuals)):
            print dihedrals[i * num_phipsi: (i + 1) * num_phipsi]
        
    for i in xrange (0, len(individuals)):
        indiv = individuals[i]
        for j in xrange (0, len(indiv.phi_dihedrals)):
            indiv.phi_dihedrals[j] = dihedrals[i * num_phipsi + j][0]
            indiv.psi_dihedrals[j] = dihedrals[i * num_phipsi + j][1]


# Generate dihedrals with BASILISK by trained DBN that interpolates continuous distribution and gives side-chain angles chis
def BasiliskSideChainDihedralsGenerator(settings, individuals):

    # MC dihedrals generation like in upper function
    # MonteCarloDihedralGenerator(settings, individuals)

    print("BASILISK_OPT BACKBONE & SIDE CHAIN DIHEDRALS ARE GENERATED")

    if (settings.verbose):
        print("ind, res, phi, psi, chi angles, ind fitnesses")

    from src.basilisk.basilisk_lib import basilisk_dbn
    from src.basilisk.basilisk_lib import basilisk_utils
    from src.basilisk.basilisk_lib import CalcAngles
    from math import degrees, radians
    
    '''
    BASILISK is a model of side chain conformational space in proteins and allows
    to sample the main degrees of freedom, the rotational chis angles, for all standard
    amino acids in continuous space. Furthermore does BASILISK incorporate a
    continuous backbone dependence, which is known to have a strong influence
    on the chis angle distributions. BASILISK has been implemented as dynamic
    Bayesian network (DBN). This package includes the fully trained DBN, an easy to use Python implementation.
    '''
    
    # Loading the residues to optimise (indexed by dihedral_residue_indexes) from the initial molecule
    mol_in = settings.initial_molecule
    res_index = []
    res_name = []
    for i in xrange(0, len(settings.dihedral_residue_indexes)):
        res = mol_in.GetResidue(settings.dihedral_residue_indexes[i])
        res_name.append(res.GetName())
        res_index.append(basilisk_utils.three_to_index(res_name[i]))
    
    # Load the DBN and compute chis and backbone angles depending on the aa
    for i in xrange(0, len(individuals)):
        indiv = individuals[i]
        
        for j in xrange(0, len(indiv.phi_dihedrals)):
            # sample a new set of angles
            # what we get back from the dbn is a tuple consisting of a
            # list of chi angles, a list with the sampled or given 
            # backbone angles and the log likelihood that was calculated 
            # for that sample 
            # with log probability in ll:
            dbn = basilisk_dbn()
            chis, bb, ll = dbn.get_sample(res_index[j], radians(indiv.phi_dihedrals[j]), radians(indiv.psi_dihedrals[j]))
            
            # without log probability
            # chis, bb, ll = dbn.get_sample(res_index[j], math.radians(indiv.phi_dihedrals[j]), math.radians(indiv.psi_dihedrals[j]), no_ll=True)
            indiv.phi_dihedrals[j] = degrees(bb[0])
            indiv.psi_dihedrals[j] = degrees(bb[1])

            chis_degrees = []
            for k in xrange(0, len(chis)):
                chis_degrees.append(degrees(chis[k]))

            indiv.chi_angles.append(chis_degrees) #list of lists containing the chi angles of each residue
            if (settings.verbose):
                print("{0}, {1}, {2}, {3}, {4}, {5}".format(i, settings.dihedral_residue_indexes[j], indiv.phi_dihedrals[j], indiv.psi_dihedrals[j], indiv.chi_angles[j], ll))

            side_chain_atoms_dict = CalcAngles.get_chi(res_name[j]) #returns a dictionary of arrays, describing the atoms for all the angles in the side chain of the current residue
            #indiv.chi_atoms.append(side_chain_atoms) #list of dictionaries
            side_chain_atoms = []
            for k in xrange(0, len(side_chain_atoms_dict)):
                side_chain_atoms.append(map(lambda x:x.upper(), side_chain_atoms_dict['x'+str(k+1)]))
                
            indiv.chi_atoms.append(side_chain_atoms) #list of 4 atoms defining each chi dihedral
