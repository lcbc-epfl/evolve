'''
The generator populates the first generation with a set of rotamers based on the selection in the input file

.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu
'''
from __future__ import print_function
from __future__ import absolute_import

# Import modules
from builtins import range
from . import Individual
from src import constants as cnts
import numpy as np

from src import MoleculeInfo as mi
from openbabel import openbabel


def initialisePopulation(settings):
    '''

    Initialisation of the first population

    Parameters
    ----------
    settings: object
        see :class:`src.JobConfig.Settings`

    Returns
    -------
    individuals: list
        a list containing objects holding each one individual with its properties and the associated molfile as OBMolec
    '''
    
    individuals = []
    
    for i in range (0, settings.population_size):
        indiv = Individual.Individual(settings)
        
        individuals.append(indiv)
        
    return individuals

def generate(settings, population, generator_op, **args):
    '''

    Generator operation that selects the coressponding generator operation function from the setting object

    Select one of the three options in the input file as follows.

    Only needed if you use composition optimization

    .. code-block:: python

        [GENERATOR]
       generator=swissprot|uniform|random


    Parameters
    ----------
    settings: object
        see :class:`src.JobConfig.Settings`
    population: list
        list of :class:`src.gaapi.Individual`

    Returns
    -------
    generator_op : function
        selected replacer

    '''
    return generator_op(settings, population, **args)

def UniformDihedralGenerator(settings, individuals):
    '''

    Generate dihedrals with uniform distribution.
    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig`
    individuals : object
        see :class:`src.gaapi.Individual`



    '''
    
    num_phipsi = len(settings.dihedral_residue_indexes)
   
    for i in range (0, settings.population_size):
        indiv = individuals[i]
        initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi, 2))
        for j in range (0, len(indiv.phi_dihedrals)):
            indiv.phi_dihedrals[j] = initial_solution[j][0]
            indiv.psi_dihedrals[j] = initial_solution[j][1]
        indiv.applyPhiPsiDihedrals(settings)


def swissprot_composition_generator(settings, individuals):
    '''

     Generates a composition distribution from the probabilities present in swissprot

    Parameters
    ----------
    settings
    individuals

    Returns
    -------

    '''
    def swissprot_to_index(min, max, selected_rotamers, selected_rotamer_types, swissprot):
        rot_count_dict = {}

        for i in range(min, max):
            rot_type = selected_rotamer_types[i]
            if rot_type in rot_count_dict.keys():
                rot_count_dict[rot_type] = rot_count_dict[rot_type] + 1
            else:
                rot_count_dict[rot_type] = 0

        rotamer_counts = rot_count_dict.values()

        rot_types = list(rot_count_dict.keys())
        # todo potentially both unbiased and swissprot should be unified
        # same code as below only that we don"t use randint but use a normalized probability
        probabilities = []
        for key in rot_types:
            probabilities.append(swissprot[key])
            pass
        # normalize probabilities to 1
        sum_probabilities = sum(probabilities)
        norm_probabilities = [float(i) / sum_probabilities for i in probabilities]
        which_type = np.random.choice(len(rot_types), p=norm_probabilities)

        chosen_rot = rot_types[which_type]

        allowed_indexes = []
        for i in range(min, max):
            rot_type = selected_rotamer_types[i]

            if (rot_type == chosen_rot):
                allowed_indexes.append(i)

        allowed_indexes = np.asarray(allowed_indexes)

        rel_index = np.random.randint(0, len(allowed_indexes))

        return allowed_indexes[rel_index]

    for i in range(0, settings.population_size):
        indiv = individuals[i]

        for j in range(0, len(settings.composition_residue_indexes)):
            indiv.composition[j] = swissprot_to_index(settings.composition_lower_bounds[j],
                                                      settings.composition_upper_bounds[j], cnts.selected_rotamers,
                                                      cnts.selected_rotamer_types, cnts.swissprot_probabilities)

        print(indiv.composition)
        indiv.applyComposition(settings)

def UniformCompositionGenerator(settings, individuals):
    '''

    This will randomly sample all the selected rotamers and thus will oversample amino acids with many rotamers (e.g Lys) compared to alanine.


    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig`
    individuals : object
        see :class:`src.gaapi.Individual`

    '''
    
    for i in range (0, settings.population_size):
        indiv = individuals[i]
        for j in range (0, len(settings.composition_residue_indexes)):
            indiv.composition[j] = np.random.randint(low=settings.composition_lower_bounds[j], high=settings.composition_upper_bounds[j])
        print(indiv.composition)
        indiv.applyComposition(settings)

        
def unbiased_protein_composition_generator(settings, individuals):
    '''

    Generates an unbiased sample across all rotamers (ie weights each rotamer type equally for initial sampling)
    All amino acids have probability of :math:`\\frac{ 1 }{n(selected amino acid types}` to be selected.

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig`
    individuals : object
        see :class:`src.gaapi.Individual`

    Returns
    -------

    '''

    def get_ubia_index(min, max, selected_rotamer_types):
        '''

        defines the indexes of allowed rotamers

        Parameters
        ----------
        min : int
            lower bound of allowed rotamers for this site
        max : int
            upper bound of allowed rotamers for this site
        selected_rotamer_types:
            list of rotamer types that are alllowed according to input file

        Returns
        -------
        allowed_indexes[rel_index]: int
            random index within bounds of this rotamer
        '''
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
            indiv.composition[j] = get_ubia_index(settings.composition_lower_bounds[j], settings.composition_upper_bounds[j], cnts.selected_rotamers)
        
        print(indiv.composition)
        indiv.applyComposition(settings)  


def MonteCarloDihedralGenerator(settings, individuals, prob_pointers=None):
    '''

    Generate dihedrals by Monte Carlo simulation

    to use this option in the inputfile use

        [MC_GENERATE]
         mc_generate_dihedrals = True
         distribution_path = 's'
         mcmove_lbound = -45.0
         mcmove_ubound = 45.0
         num_mc_steps = 1000
         dihedral_probability_pointers = 1 2 3 1

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig`
    individuals : object
        see :class:`src.gaapi.Individual`
    prob_pointers

    Returns
    -------

    '''

    
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



def BasiliskSideChainDihedralsGenerator(settings, individuals):
    '''

    Generate dihedrals with BASILISK by trained DBN that interpolates continuous distribution and gives side-chain angles chis
    MC dihedrals generation like in upper function
    MonteCarloDihedralGenerator(settings, individuals)

    Parameters
    ----------
    settings
    individuals

    Returns
    -------

    '''


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
    '''

    #TODO Nick

    Parameters
    ----------
    obmol

    Returns
    -------

    '''
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
    '''

    #TODO Nick

    Parameters
    ----------
    obmol

    Returns
    -------

    '''
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

    # TODO Nick
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
    '''
     Test for standalone run for the basilisk test case. 
     Should be moved to a test #TODO
    '''
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
    
