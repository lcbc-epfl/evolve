'''
Individual sets up the Invidiuals for EVOLVE computation

Depending on context the individual is deepcopied or initalised.

.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu
'''
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
from src import MoleculeInfo as mi
from src import MoleculeCreator as mcr
import copy
from openbabel import openbabel
import numpy as np
import pygmo as pg

class Individual(object):
    
    def __init__(self, settings, orig=None):
        if orig is None:
            self.initialise_constructor(settings)
        else:
            self.copy_constructor(settings, orig)
        
    def initialise_constructor(self, settings):
        '''

        Initalises the Individual from the read in PDB file and initializes all variables.

        Parameters
        ----------
        settings : object
            see :class:`src.JobConfig.Settings`

        Returns
        -------

        '''
        
        if (settings.multi_individual):
            self.mol=[None]*settings.no_frames
            for i,x in enumerate(settings.initial_molecule):
                self.mol[i]=openbabel.OBMol(x)
                pass
            pass
        else:
            self.mol = openbabel.OBMol(settings.initial_molecule)
            pass

        self.stability = 0
        # print self.mol, settings.initial_molecule
        self.dihedral_residue_indexes = None
        self.phi_dihedrals = None
        self.psi_dihedrals = None
        
        self.composition_residue_indexes = None
        self.composition = None

        self.chi_angles = None
        self.chi_atoms = None
        
        self.energies = {"wt":None, "wt_corr":None, "mut":None, "mut_corr":None}
        
        
        if (settings.backbone_dihedral_optimization):
            
            self.dihedral_residue_indexes = settings.dihedral_residue_indexes
            
            self.phi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
           
            self.updatePhiPsiDihedrals(settings)
            
        if (settings.composition_optimization):
            
            self.composition_residue_indexes = settings.composition_residue_indexes
            self.composition = np.zeros(len(self.composition_residue_indexes), dtype=np.int)
        
        if (settings.sidechain_dihedral_optimisation):
            self.updateChiAngles(settings)

            
        self.fitnesses = np.zeros(len(settings.evaluators))
    
    # doing this to get around PySwigObject's  inability to be deepcopied  
    def copy_constructor(self, settings, orig):
        import copy

        if (settings.multi_individual):
            self.mol=[None]*settings.no_frames
            for i,x in enumerate(orig.mol):
                self.mol[i]=openbabel.OBMol(x)
                pass
            pass
        else:
            self.mol = openbabel.OBMol(orig.mol)
            pass
        
        # print self.mol, orig.mol
        self.dihedral_residue_indexes = copy.deepcopy(orig.dihedral_residue_indexes)
        
        if (self.dihedral_residue_indexes != None):
            self.phi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            
            for i in range (0, len(self.dihedral_residue_indexes)):
                self.phi_dihedrals[i] = orig.phi_dihedrals[i]
                self.psi_dihedrals[i] = orig.psi_dihedrals[i]
        
        
        self.composition_residue_indexes = copy.deepcopy(orig.composition_residue_indexes)
        
        if (self.composition_residue_indexes != None):
            self.composition = np.zeros(len(self.composition_residue_indexes), dtype=np.int)
            for i in range (0, len(self.composition_residue_indexes)):
                self.composition[i] = orig.composition[i]

        if (settings.sidechain_dihedral_optimisation):
            self.chi_angles = copy.deepcopy(orig.chi_angles)
        
        
        self.fitnesses = copy.deepcopy(orig.fitnesses)
        self.energies = orig.energies
    
    def saveMol(self, file_path):
        '''
        saves molecule as file specified by the path that was passed

        Parameters
        ----------
        file_path: str
            path to the file including ending to retrieve format. Must be compatible to openbabel output

        '''
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat(file_path.split('.')[1])
        
        obconv.WriteFile(self.mol, file_path)
        
    def collectChiAtomPointers(self):
        '''

        #TODO
        Nick

        Returns
        -------

        '''
        chi_atoms = []
            
         # OBAtom references are now different -need to refind corresponding atoms
        for i in range (self.mol.NumResidues()):
            
            obres = self.mol.GetResidue(i)
            sidechain_atoms_dict = mi.get_chi_atoms(obres)
        
            if (not sidechain_atoms_dict or obres.GetName() == "PRO"):
                chi_atoms.append([])
                continue
            
            num_chi_angles = len(sidechain_atoms_dict)
            # print num_chi_angles
            res_chi_atoms = []
            
            for j in range(num_chi_angles):
                atoms = [None, None, None, None]
        
                sidechain_torsion_atomsID = [atomID.upper() for atomID in sidechain_atoms_dict['x' + str(j + 1)]]
                
                (IDatoms_in_res, OBatoms_in_res, NUMatoms_in_res) = mi.get_atoms(obres)
                
                for k in range(0, 4):
                    atoms[k] = OBatoms_in_res[IDatoms_in_res.index(sidechain_torsion_atomsID[k])]
        
        
                res_chi_atoms.append(atoms)
                
            chi_atoms.append(res_chi_atoms)
            
        return chi_atoms

    def setStability(self,value):
        '''

        Set stability to defined values

        Parameters
        ----------
        value

        Returns
        -------

        '''
        
        self.stability=value

    def setEnergies(self,wt, wt_corr, mut, mut_corr):
        '''

        Set energies 

        Parameters
        ----------
        value

        Returns
        -------

        '''
        self.energies.update({"wt":wt, "wt_corr":wt_corr, "mut":mut, "mut_corr":mut_corr})

    def updatePhiPsiDihedrals(self, settings):
        '''

        #TODO
        Nick

        Parameters
        ----------
        settings

        Returns
        -------

        '''
        phi_psi_dihedrals = self.getPhiPsiDihedrals(settings)
        
        for j in range (0, len(self.dihedral_residue_indexes)):
            self.phi_dihedrals[j] = phi_psi_dihedrals[j][0]
            self.psi_dihedrals[j] = phi_psi_dihedrals[j][1]
                
    def getPhiPsiDihedrals(self, settings):
        '''

        #TODO Nick

        Parameters
        ----------
        settings

        Returns
        -------

        '''
        '''returns (phi, psi) tuple'''
        if (self.dihedral_residue_indexes == None):
            return None
               
        return mi.getPhiPsiDihedralByAtomIndex(self.mol, settings.backbone_dihedral_atom_idxs)

    def applyPhiPsiDihedrals(self, settings):
        '''
        #TODO Nick
        '''
        if (self.phi_dihedrals is None):
            return 
        
        if (settings.backbone_dihedral_atom_idxs is None or len(settings.backbone_dihedral_atom_idxs) == 0):
            return

        for j in range(len(settings.dihedral_residue_indexes)):
            
            # print j, self.mol.GetResidue(j).GetName()
            
            if (len(settings.backbone_dihedral_atom_idxs[j][0]) == 0):
                # no phi angle
                continue
            
            # set phi
            phi_angle_atoms = settings.backbone_dihedral_atom_idxs[j][0]
            self.mol.SetTorsion(phi_angle_atoms[0], phi_angle_atoms[1], phi_angle_atoms[2], phi_angle_atoms[3], np.radians(self.phi_dihedrals[j]))
            
            if (len(settings.backbone_dihedral_atom_idxs[j][1]) == 0):
                # no psi angle
                continue
            
            # set psi
            psi_angle_atoms = settings.backbone_dihedral_atom_idxs[j][1]
            self.mol.SetTorsion(psi_angle_atoms[0], psi_angle_atoms[1], psi_angle_atoms[2], psi_angle_atoms[3], np.radians(self.psi_dihedrals[j]))

            
            
    def applyComposition(self, settings):
        '''

        selects rotamer for all positions to mutate based on input file
        for this takes the mol2 file of the residue to mutate and uses the swapsidechain  method to exchange sidechains

        run for every iteration of the ga after mutate and crossover publications

        Parameters
        ----------
        settings:



        '''
     
        from src import constants

        if (self.composition_residue_indexes == None):
            return
        print(self.composition_residue_indexes)
        cmp = []
        if (settings.verbose):
            print("Created Individual:")
        for i in range(0, len(self.composition_residue_indexes)):
            if (settings.multi_individual):
                for x in range(0, settings.no_frames):
                    # Need to include frag inside the loop, because the open babel objects get modified inside swapsidechain
                    frag = openbabel.OBMol()
                    obConversion = openbabel.OBConversion()
                    obConversion.SetInFormat("mol2")

                    rotamer_type = constants.rotamers[self.composition[i]]

                    obConversion.ReadFile(frag, settings.composition_library + "/" + rotamer_type + ".mol2")
                    try:
                        mcr.swapsidechain(settings, self.mol[x], self.composition_residue_indexes[i], frag)
                    except AttributeError as e:
                        print(f"{i} went wrong, ignore it")
                    pass
                cmp.append(rotamer_type) # append only after loop because all individuals are identical
            else:
                frag = openbabel.OBMol()
                obConversion = openbabel.OBConversion()
                obConversion.SetInFormat("mol2")

                rotamer_type = constants.selected_rotamers[self.composition[i]]
                cmp.append(rotamer_type)
                if (settings.verbose):
                    print(cmp)
                obConversion.ReadFile(frag, settings.composition_library + "/" + rotamer_type + ".mol2")
                mcr.swapsidechain(settings, self.mol, self.composition_residue_indexes[i], frag)


    


    def applyChiDihedrals(self, settings):
        '''

        set dihedrals for dihedral optimization

        Parameters
        ----------
        settings

        Returns
        -------

        '''

        if (self.chi_angles is None or len(self.chi_angles) == 0):
            return

        for j in range(len(settings.chi_dihedral_atom_idxs)):
            
            # print j, self.mol.GetResidue(j).GetName()
            
            if (len(settings.chi_dihedral_atom_idxs[j]) == 0):
                # residue has no chi atoms, or needs to be ignored (PRO)
                continue

            for k in range(0, len(settings.chi_dihedral_atom_idxs[j])):
                
                chi_angle_atoms = settings.chi_dihedral_atom_idxs[j][k]
        
                self.mol.SetTorsion(chi_angle_atoms[0], chi_angle_atoms[1], chi_angle_atoms[2], chi_angle_atoms[3], np.radians(self.chi_angles[j][k]))
                
    def updateChiAngles(self, settings):
       self.chi_angles = mi.getChiDihedralsByAtomIndex(self.mol, settings.chi_dihedral_atom_idxs)

     
    def dominates(self, individual):
        '''

        to compare individuals. For one fitness the minimal value is used.
        For multiple fitnesses non_domination rank is used as implemented in the pygmo library.

        Parameters
        ----------
        individual : object
            see :class:`src.gaapi.Individual`

        Returns
        -------
        domination : bool
            individual 1 dominates individual 2
        '''
        if (len(self.fitnesses) > 1) :
            print("Dominance:", self.fitnesses, individual.fitnesses, pg.pareto_dominance(obj1=self.fitnesses, obj2=individual.fitnesses))
            return pg.pareto_dominance(obj1 = self.fitnesses,obj2=individual.fitnesses)
            pass
        else:
            return (self.fitnesses[0] < individual.fitnesses[0])
        
    def setFitness(self, index, fitness):
        self.fitnesses[index] = fitness
        
    def getFitness(self, index):
        return self.fitnesses[index]
        
    def __lt__(self, other):
        '''
        use dominate function for comparisons as the individual can have multiple fitnesses
        Parameters
        ----------
        other: object
            see :class:`src.gaapi.Individual`


        Returns
        -------
        dominiation: bool
            1 if self < other is true, 0 if self <0 false

        '''
        return self.dominates(other)
    
    def __gt__(self, other):
        '''
            use dominate function for comparisons as the individual can have multiple fitnesses
            Parameters
            ----------
            other: object
                see :class:`src.gaapi.Individual`


            Returns
            -------
            dominiation: bool
                1 if self > other is true, 0 if self >0 false

            '''
        return other.dominates(self)
    
class PSOIndividual(Individual):
    '''

    #TODO
    NIck, perhaps can be deleted


    '''
    def __init__(self, settings, orig=None):
        if orig is None:
            self.initialise_constructor(settings)
            
            self.phi_velocity_vector = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_velocity_vector = np.zeros(len(self.dihedral_residue_indexes))
            
            self.topological_neighbours = None
            self.topological_best_neighbour = None
            
            self.personal_best_phi = np.zeros(len(self.dihedral_residue_indexes))
            self.personal_best_psi = np.zeros(len(self.dihedral_residue_indexes))
            self.personal_best_fitness = None
            
        else:
            self.copy_constructor(settings, orig)
            self.phi_velocity_vector = orig.phi_velocity_vector
            self.psi_velocity_vector = orig.psi_velocity_vector
            
            self.topological_neighbours = orig.topological_neighbours
            self.topological_best_neighbour = orig.topological_best_neighbour
            
            self.personal_best_phi = orig.personal_best_phi
            self.personal_best_psi = orig.personal_best_psi
            
            
    def updateVelocity(self, population, weight, c1, c2):
        
        r1 = np.random.rand((len(self.phi_dihedrals)))
        r2 = np.random.rand((len(self.phi_dihedrals)))
        
        self.phi_velocity_vector = (weight * self.phi_velocity_vector) + (c1 * np.multiply(self.personal_best_phi - self.phi_dihedrals, r1)) + \
        (c2 * np.multiply(population[self.topological_best_neighbour].phi_dihedrals - self.phi_dihedrals, r2))
        
        r1 = np.random.rand((len(self.psi_dihedrals)))
        r2 = np.random.rand((len(self.psi_dihedrals)))

        self.psi_velocity_vector = (weight * self.psi_velocity_vector) + (c1 * np.multiply(self.personal_best_psi - self.psi_dihedrals, r1)) + \
        (c2 * np.multiply(population[self.topological_best_neighbour].psi_dihedrals - self.psi_dihedrals, r2))
        
        for i in range(0, len(self.phi_dihedrals)):
            
            if (self.phi_velocity_vector[i] < -45):
                self.phi_velocity_vector[i] = -45
            if (self.phi_velocity_vector[i] > 45):
                self.phi_velocity_vector[i] = 45
                
            if (self.psi_velocity_vector[i] < -45):
                self.psi_velocity_vector[i] = -45
            if (self.psi_velocity_vector[i] > 45):
                self.psi_velocity_vector[i] = 45
    
    
    def updatePosition(self):
        self.phi_dihedrals = self.phi_dihedrals + self.phi_velocity_vector
        self.psi_dihedrals = self.psi_dihedrals + self.psi_velocity_vector
        
        for i in range(0, len(self.phi_dihedrals)):
            if (self.phi_dihedrals[i] < -180):
                self.phi_dihedrals[i] = self.phi_dihedrals[i] + 360
            if (self.phi_dihedrals[i] > 180):
                self.phi_dihedrals[i] = self.phi_dihedrals[i] - 360
                
            if (self.psi_dihedrals[i] < -180):
                self.psi_dihedrals[i] = self.psi_dihedrals[i] + 360
            if (self.psi_dihedrals[i] > 180):
                self.psi_dihedrals[i] = self.psi_dihedrals[i] - 360
        
    def updatePersonalBest(self):
        if self.fitnesses[0] < self.personal_best_fitness:
            self.personal_best_fitness = self.fitnesses[0]
            self.personal_best_phi = phi_dihedrals
            self.personal_best_psi = psi_dihedrals
        
    
    def updateTopologicalBestNeighbour(self, population):
        for i in range(len(self.topological_neighbours)):
            if (population[self.topological_neighbours[i]].dominates(population[self.topological_best_neighbour])):
                self.topological_best_neighbour = self.topological_neighbours[i]
    
    
