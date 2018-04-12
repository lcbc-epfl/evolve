'''
main.py

@author: Nicholas Browning
'''
from src import MoleculeInfo as mi
from src import MoleculeCreator as mcr
import copy
import openbabel
import numpy as np

class Individual(object):
    
    def __init__(self, settings, orig=None):
        self.init = False  # True if indiv is the initial structure of GA
        if orig is None:
            self.initialise_constructor(settings)
        else:
            self.copy_constructor(settings, orig)
        
    def initialise_constructor(self, settings):
        
        self.mol = openbabel.OBMol(settings.initial_molecule)
        
        self.dihedral_residue_indexes = None
        self.phi_dihedrals = None
        self.psi_dihedrals = None
        
        self.composition_residue_indexes = None
        self.composition = None

        self.chi_angles = None
        self.chi_atoms = None
        
        
        if (settings.dihedral_optimization):
            
            self.dihedral_residue_indexes = settings.dihedral_residue_indexes
            
            self.phi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
           
            phi_psi_dihedrals = self.getPhiPsiDihedrals()
           
            for i in xrange (0, len(phi_psi_dihedrals)):
                self.phi_dihedrals[i] = phi_psi_dihedrals[i][0]
                self.psi_dihedrals[i] = phi_psi_dihedrals[i][1]
            
        if (settings.composition_optimization):
            
            self.composition_residue_indexes = settings.composition_residue_indexes
            self.composition = np.zeros(len(self.composition_residue_indexes), dtype=np.int)

        if (settings.basilisk_and_sidechains):

            self.chi_angles = []  # Sidechain chi dihedrals
            self.chi_atoms = []  # Sidechain atoms defining chi dihedrals, name not OBAtom

            
        self.fitnesses = np.zeros(len(settings.evaluators))
    
    # doing this to get around PySwigObject's  inability to be deepcopied  
    def copy_constructor(self, settings, orig):
        
        self.mol = openbabel.OBMol(orig.mol)
        
        self.dihedral_residue_indexes = orig.dihedral_residue_indexes
        
        if (self.dihedral_residue_indexes != None):
            self.phi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            
            for i in xrange (0, len(self.dihedral_residue_indexes)):
                self.phi_dihedrals[i] = orig.phi_dihedrals[i]
                self.psi_dihedrals[i] = orig.psi_dihedrals[i]
        
        self.composition_residue_indexes = orig.composition_residue_indexes
        
        if (self.composition_residue_indexes != None):
            self.composition = np.zeros(len(self.composition_residue_indexes), dtype=np.int)
            for i in xrange (0, len(self.composition_residue_indexes)):
                self.composition[i] = orig.composition[i]

        self.chi_angles = orig.chi_angles
        
        if (self.chi_angles != None):
            self.chi_angles = orig.chi_angles
            self.chi_atoms = orig.chi_atoms
            
        
        self.fitnesses = np.zeros(len(orig.fitnesses))
        
        for i in xrange (0, len(orig.fitnesses)):
            self.fitnesses[i] = orig.fitnesses[i]
        

    def getPhiPsiDihedrals(self):
        '''returns (phi, psi) tuple'''
        if (self.dihedral_residue_indexes == None):
            return None
        
        return mi.getPhiPsiDihedrals(self.mol, self.dihedral_residue_indexes)

    def applyPhiPsiDihedrals(self):
        if (self.phi_dihedrals is None):
            return
        
        for i in xrange (0, len(self.dihedral_residue_indexes)):
            obres = self.mol.GetResidue(self.dihedral_residue_indexes[i])
            mi.SetPhiPsi(self.mol, obres, self.phi_dihedrals[i], self.psi_dihedrals[i])
            
            
    def applyComposition(self, settings):
        from src import constants
        
        if (self.composition_residue_indexes == None):
            return

        for i in xrange(0, len(self.composition_residue_indexes)):
            frag = openbabel.OBMol()
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("mol2")
            
            rotamer_type = constants.rotamers[self.composition[i]]
            
            obConversion.ReadFile(frag, settings.composition_library + "/" + rotamer_type + ".mol2")
            
            # print self.composition_residue_indexes[i], rotamer_type + ".mol2"
            mcr.swapsidechain(self.mol, self.composition_residue_indexes[i], frag)

    
    def get_chi_dihedrals_per_res(self):

        if (self.chi_angles == None or self.chi_angles == []):
            return
        
        (sorted_atoms_ID, sorted_atoms_OB) = mi.get_atoms_per_residue(self.mol)

        ss_angles = [[] for _ in xrange(0, len(self.dihedral_residue_indexes))]

        for j in self.dihedral_residue_indexes:

            for k in xrange(0, len(self.chi_angles[j])):
                torsion_atom_ID = self.chi_atoms[j][k]
                # print(torsion_atom_ID)
                # print(sorted_atoms_ID[j])
                
                torsion_atom_locs = []
                torsion_ob_atoms = []

                for l in xrange(0, 4):
                    torsion_atom_locs.append(sorted_atoms_ID[j].index(torsion_atom_ID[l]))
                    torsion_ob_atoms.append(sorted_atoms_OB[j][torsion_atom_locs[l]])

                ss_angles[j].append(self.mol.GetTorsion(torsion_ob_atoms[0], torsion_ob_atoms[1], torsion_ob_atoms[2], torsion_ob_atoms[3]))

        return ss_angles
         

    def apply_chi_dihedrals(self):

        if (self.chi_angles is None or len(self.chi_angles) == 0):
            return
        
        (sorted_atoms_ID, sorted_atoms_OB) = mi.get_atoms_per_residue(self.mol)

        for j in self.dihedral_residue_indexes:

            for k in xrange(0, len(self.chi_angles[j])):
                torsion_atom_ID = self.chi_atoms[j][k]
                # print(torsion_atom_ID)
                # print(sorted_atoms_ID[j])
                
                torsion_atom_locs = []
                torsion_ob_atoms = []

                for l in xrange(0, 4):
                    torsion_atom_locs.append(sorted_atoms_ID[j].index(torsion_atom_ID[l]))
                    torsion_ob_atoms.append(sorted_atoms_OB[j][torsion_atom_locs[l]])

                self.mol.SetTorsion(torsion_ob_atoms[0], torsion_ob_atoms[1], torsion_ob_atoms[2], torsion_ob_atoms[3], np.radians(self.chi_angles[j][k]))
                
                
    def dominates(self, individual):
        if (len(self.fitnesses) > 1) :
            # TODO
            pass
        else:
            return (self.fitnesses[0] < individual.fitnesses[0])
        
    def setFitness(self, index, fitness):
        self.fitnesses[index] = fitness
        
    def __lt__(self, other):
        return self.dominates(other)
    
    def __gt__(self, other):
        return other.dominates(self)
    
class PSOIndividual(Individual):
    
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
    
    
