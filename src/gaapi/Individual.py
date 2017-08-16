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
        
        self.dihedral_residue_indexes = None
        self.phi_dihedrals = None
        self.psi_dihedrals = None
        
        self.composition_residue_indexes = None
        self.composition = None
        
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
        
        self.fitnesses = np.zeros(len(orig.fitnesses))
        
        for i in xrange (0, len(orig.fitnesses)):
            self.fitnesses[i] = orig.fitnesses[i]
        

    def getPhiPsiDihedrals(self):
        '''returns (phi, psi) tuple'''
        if (self.dihedral_residue_indexes == None):
            return None
        
        return mi.getPhiPsiDihedrals(self.mol, self.dihedral_residue_indexes)

    def applyPhiPsiDihedrals(self, settings):
        if (self.phi_dihedrals == None):
            return
        
        for i in xrange (0, len(self.dihedral_residue_indexes)):
            obres = self.mol.GetResidue(self.dihedral_residue_indexes[i])
            mi.SetPhiPsi(self.mol, obres, self.phi_dihedrals[i], self.psi_dihedrals[i])
            
    def applyComposition(self, settings):
        from src import constants
        
        if (self.composition_residue_indexes == None):
            return

        for i in xrange (0, len(self.composition_residue_indexes)):
            frag = openbabel.OBMol()
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("mol2")
            
            rotamer_type = constants.rotamers[self.composition[i]]
            
            obConversion.ReadFile(frag, settings.composition_library + "/" + rotamer_type + ".mol2") 
            
            # print self.composition_residue_indexes[i], rotamer_type + ".mol2"
            mcr.swapsidechain(self.mol, self.composition_residue_indexes[i], frag)
            
        
            
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
            
    
    
    
