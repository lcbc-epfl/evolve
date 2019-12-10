'''
generators.py

@author: Nicholas Browning
'''
from __future__ import print_function
from __future__ import absolute_import

# Import modules
from builtins import range
from . import Individual
from src import constants as cnts
import numpy as np

from src import MoleculeInfo as mi
import openbabel


# Initialisation of the first population
def initialisePopulation(settings):
    
    individuals = []
    
    for i in range (0, settings.population_size):
        indiv = Individual.Individual(settings)
        
        individuals.append(indiv)
        
    return individuals


# Generate dihedrals with uniform distribution
def UniformDihedralGenerator(settings, individuals):
    
    num_phipsi = len(settings.dihedral_residue_indexes)
   
    for i in range (0, settings.population_size):
        indiv = individuals[i]
        initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi, 2))
        for j in range (0, len(indiv.phi_dihedrals)):
            indiv.phi_dihedrals[j] = initial_solution[j][0]
            indiv.psi_dihedrals[j] = initial_solution[j][1]
        indiv.applyPhiPsiDihedrals(settings)
    

def UniformCompositionGenerator(settings, individuals):
    
    for i in range (0, settings.population_size):
        indiv = individuals[i]
        for j in range (0, len(settings.composition_residue_indexes)):
            indiv.composition[j] = np.random.randint(low=settings.composition_lower_bounds[j], high=settings.composition_upper_bounds[j])
        print(indiv.composition)
        indiv.applyComposition(settings)

        
def unbiased_protein_composition_generator(settings, individuals):
    
    # Generates an unbiased sample across all rotamers (ie weights each rotamer type equally for initial sampling)
    def get_ubia_index(min, max, selected_rotamers, selected_rotamer_types):
        rot_count_dict = {}

        for i in range(min, max):
            rot_type = selected_rotamer_types[i]
            if rot_type in list(rot_count_dict.keys()):
                rot_count_dict[rot_type] = rot_count_dict[rot_type] + 1
            else:
                rot_count_dict[rot_type] = 0
    
        rotamer_counts = list(rot_count_dict.values())
        rot_types = list(rot_count_dict.keys())
        
        which_type = np.random.randint(low=0, high=len(rot_types))
        
        chosen_rot = rot_types[which_type]
        
        allowed_indexes = []
        
        for i in range(min, max):
            rot_type = selected_rotamer_types[i]

            if (rot_type == chosen_rot):
                allowed_indexes.append(i)
        
        allowed_indexes = np.asarray(allowed_indexes)
                
        rel_index = np.random.randint(0, len(allowed_indexes))
        
        return allowed_indexes[rel_index]
    
    for i in range (0, settings.population_size):
        indiv = individuals[i]
        
        for j in range (0, len(settings.composition_residue_indexes)):
            indiv.composition[j] = get_ubia_index(settings.composition_lower_bounds[j], settings.composition_upper_bounds[j], cnts.selected_rotamers, cnts.selected_rotamer_types)
        
        print(indiv.composition)
        indiv.applyComposition(settings)  

            
# Generate dihedrals by Monte Carlo
def MonteCarloDihedralGenerator(settings, individuals, prob_pointers=None):
    
#     [MC_GENERATE]
# mc_generate_dihedrals = True
# distribution_path = 's'
# mcmove_lbound = -45.0
# mcmove_ubound = 45.0
# num_mc_steps = 1000
# dihedral_probability_pointers = 1 2 3 1
    
    from deprecated.montecarlo import dihedral_opt as mcopt
 
    probs = mcopt.loadDefaultProbabilityDistributions()
    
    print(len(probs))
    num_phipsi = len(settings.dihedral_residue_indexes)
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi * len(individuals), 2))
    
    print("INITIAL DIHEDRALS UNIFORMLY AND RANDOMLY GENERATED")
            
    move_lows = [-30] * num_phipsi * settings.population_size
    move_highs = [30] * num_phipsi * settings.population_size
    
    prob_idxs = []
    
    if prob_pointers == None:
        for i in range(0, len(individuals)):
            for v in getMCProbPinters(individuals[i].mol): 
                prob_idxs.append(v)
    else:
        prob_idxs = prob_pointers
    
    dihedrals = mcopt.optimisePhiPsi(probs, prob_idxs, initial_solution, 1000, move_lows, move_highs)
    print("MC_OPTIMIZED DIHEDRALS GENERATED")
        
    for i in range (0, len(individuals)):
        indiv = individuals[i]
        for j in range (0, len(indiv.phi_dihedrals)):
            indiv.phi_dihedrals[j] = dihedrals[i * num_phipsi + j][0]
            indiv.psi_dihedrals[j] = dihedrals[i * num_phipsi + j][1]


# Generate dihedrals with BASILISK by trained DBN that interpolates continuous distribution and gives side-chain angles chis
def BasiliskSideChainDihedralsGenerator(settings, individuals):

    # MC dihedrals generation like in upper function
    # MonteCarloDihedralGenerator(settings, individuals)

    print("BASILISK_OPT BACKBONE & SIDE CHAIN DIHEDRALS ARE GENERATED")

    if (settings.verbose):
        print("ind, res, phi, psi, chi angles, ind fitnesses")

    from deprecated.basilisk.basilisk_lib import basilisk_dbn
    from deprecated.basilisk.basilisk_lib import basilisk_utils
    from deprecated.basilisk.basilisk_lib import CalcAngles
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
    for i in range(0, len(settings.dihedral_residue_indexes)):
        res = mol_in.GetResidue(settings.dihedral_residue_indexes[i])
        res_name.append(res.GetName())
        res_index.append(basilisk_utils.three_to_index(res_name[i]))
        
    # print (res_index)
    # Load the DBN and compute chis and backbone angles depending on the aa
    for i in range(0, len(individuals)):
        indiv = individuals[i]
        
        for j in range(0, len(indiv.phi_dihedrals)):
           
            if len(indiv.chi_angles[j]) == 0:
               continue
            
            dbn = basilisk_dbn()
 
            chis, bb, ll = dbn.get_sample(res_index[j])  # Nick fix - call above doesn't provide independant BB samples!
            
            indiv.phi_dihedrals[j] = degrees(bb[0]) - 180.0  # range is from 0 - 360 by default
            indiv.psi_dihedrals[j] = degrees(bb[1]) - 180.0

            for k in range(0, len(chis)):
                indiv.chi_angles[j][k] = degrees(chis[k]) - 180.0

        
def getBasiliskSample(obmol):
    from deprecated.basilisk.basilisk_lib import basilisk_dbn
    
    res_name, res_index = getResidueIndexes(obmol)
    
    dbn = basilisk_dbn()
    
    chis = []
    bbs = []
    lls = []
    
    for residx in res_index:
        chi, bb, ll = dbn.get_sample(residx)
        
        chis.append(chi)
        bbs.append(bb)
        lls.append(ll)
    
    return chis, bbs, lls
    
       
def getResidueIndexes(obmol):
    from deprecated.basilisk.basilisk_lib import basilisk_utils
    from src import MoleculeInfo as mi
    res_index = []
    res_name = []
    for i in range(0, obmol.NumResidues()):
        res = obmol.GetResidue(i)
        res_name.append(res.GetName())
        res_index.append(basilisk_utils.three_to_index(res.GetName()))
        
    return res_name, res_index


def getMCProbPinters(obmol):
    
    pointers = []
    for i in range(obmol.NumResidues()):
        
        res = obmol.GetResidue(i)
        
        resType = res.GetName()
        
        if (resType == "PRO"):
            pointers.append(3)
        elif (i + 1 < obmol.NumResidues() and obmol.GetResidue(i + 1).GetName() == "PRO"):
            pointers.append(2)
        elif (resType == "GLY"):
            pointers.append(1)
        else:
            pointers.append(0)
            
    return pointers


if __name__ == "__main__":
    import openbabel
    import sys
    from math import degrees, radians
    from deprecated.montecarlo import dihedral_opt as mcopt
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, sys.argv[1])
    
    print(getResidueIndexes(mol))
    chis, bbs, lls = getBasiliskSample(mol)
    
    for i, v in enumerate(bbs):
        obres = mol.GetResidue(i)
        mi.SetPhiPsi(mol, obres, degrees(v[0]), degrees(v[1]))
    
    obConversion.WriteFile(mol, "test1.pdb")
    
    from src.montecarlo import dihedral_opt as mcopt
 
    probs = mcopt.loadDefaultProbabilityDistributions()
    
    num_phipsi = 20
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi , 2))
            
    move_lows = [-45] * num_phipsi
    move_highs = [-45] * num_phipsi
    
    dihedrals = mcopt.optimisePhiPsi(probs, getMCProbPinters(mol), initial_solution, 2000, move_lows, move_highs)
    
    for i, v in enumerate(bbs):
        obres = mol.GetResidue(i)
        mi.SetPhiPsi(mol, obres, dihedrals[i][0], dihedrals[i][1])
    
    obConversion.WriteFile(mol, "test2.pdb")
    
    print(dihedrals)
    
