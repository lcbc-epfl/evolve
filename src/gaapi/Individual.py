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
        if orig is None:
            self.initialise_constructor(settings)
        else:
            self.copy_constructor(settings, orig)
        
    def initialise_constructor(self, settings):
        
        self.mol = openbabel.OBMol(settings.initial_molecule)
        
        # print self.mol, settings.initial_molecule
        self.dihedral_residue_indexes = None
        self.phi_dihedrals = None
        self.psi_dihedrals = None
        
        self.composition_residue_indexes = None
        self.composition = None

        self.chi_angles = None
        self.chi_atoms = None
        
        
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
        self.mol = openbabel.OBMol(orig.mol)
        
        # print self.mol, orig.mol
        self.dihedral_residue_indexes = copy.deepcopy(orig.dihedral_residue_indexes)
        
        if (self.dihedral_residue_indexes != None):
            self.phi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            self.psi_dihedrals = np.zeros(len(self.dihedral_residue_indexes))
            
            for i in xrange (0, len(self.dihedral_residue_indexes)):
                self.phi_dihedrals[i] = orig.phi_dihedrals[i]
                self.psi_dihedrals[i] = orig.psi_dihedrals[i]
        
        
        self.composition_residue_indexes = copy.deepcopy(orig.composition_residue_indexes)
        
        if (self.composition_residue_indexes != None):
            self.composition = np.zeros(len(self.composition_residue_indexes), dtype=np.int)
            for i in xrange (0, len(self.composition_residue_indexes)):
                self.composition[i] = orig.composition[i]

        if (settings.sidechain_dihedral_optimisation):
            self.chi_angles = copy.deepcopy(orig.chi_angles)
        
        
        self.fitnesses = copy.deepcopy(orig.fitnesses)
    
    def saveMol(self, file_path):
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat(file_path.split('.')[1])
        
        obconv.WriteFile(self.mol, file_path)
        
    def collectChiAtomPointers(self):
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
        
    def updatePhiPsiDihedrals(self, settings):
        phi_psi_dihedrals = self.getPhiPsiDihedrals(settings)
        
        for j in xrange (0, len(self.dihedral_residue_indexes)):
            self.phi_dihedrals[j] = phi_psi_dihedrals[j][0]
            self.psi_dihedrals[j] = phi_psi_dihedrals[j][1]
                
    def getPhiPsiDihedrals(self, settings):
        '''returns (phi, psi) tuple'''
        if (self.dihedral_residue_indexes == None):
            return None
               
        return mi.getPhiPsiDihedralByAtomIndex(self.mol, settings.backbone_dihedral_atom_idxs)

    def applyPhiPsiDihedrals(self, settings):
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
     
        from src import constants
        
        if (self.composition_residue_indexes == None):
            return

        for i in xrange(0, len(self.composition_residue_indexes)):
            frag = openbabel.OBMol()
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("mol2")
            
            rotamer_type = constants.rotamers[self.composition[i]]
            
            obConversion.ReadFile(frag, settings.composition_library + "/" + rotamer_type + ".mol2")
        
            mcr.swapsidechain(self.mol, self.composition_residue_indexes[i], frag)


    def applyChiDihedrals(self, settings):

        if (self.chi_angles is None or len(self.chi_angles) == 0):
            return

        for j in range(len(settings.chi_dihedral_atom_idxs)):
            
            # print j, self.mol.GetResidue(j).GetName()
            
            if (len(settings.chi_dihedral_atom_idxs[j]) == 0):
                # residue has no chi atoms, or needs to be ignored (PRO)
                continue

            for k in xrange(0, len(settings.chi_dihedral_atom_idxs[j])):
                
                chi_angle_atoms = settings.chi_dihedral_atom_idxs[j][k]
        
                self.mol.SetTorsion(chi_angle_atoms[0], chi_angle_atoms[1], chi_angle_atoms[2], chi_angle_atoms[3], np.radians(self.chi_angles[j][k]))
                
    def updateChiAngles(self, settings):
       self.chi_angles = mi.getChiDihedralsByAtomIndex(self.mol, settings.chi_dihedral_atom_idxs)

     
    def dominates(self, individual):
        if (len(self.fitnesses) > 1) :
            # TODO
            pass
        else:
            return (self.fitnesses[0] < individual.fitnesses[0])
        
    def setFitness(self, index, fitness):
        self.fitnesses[index] = fitness
        
    def getFitness(self, index):
        return self.fitnesses[index]
        
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
    
    
