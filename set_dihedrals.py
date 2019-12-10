'''
printMoleculeInfo

@author: Nicholas Browning
'''
from __future__ import print_function
# Use this to print residue information
from builtins import range
import argparse
from src import MoleculeInfo as mi
import openbabel
import numpy as np


def printResidueInfo(mol):
    print("Num Atoms:", mol.NumAtoms())
    print("Num Bonds:", mol.NumBonds())
    print("Num Residues:", mol.NumResidues())
    
    
    print("Phi, Psi Dihedrals:")
    dihedrals = mi.getAllPhiPsiDihedrals(mol)
    
    for i in range (0, mol.NumResidues()):
        res = mol.GetResidue(i)
        
        print("Res", i, res.GetName())
        
        print("[phi, psi]", dihedrals[i][0], dihedrals[i][1])
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    parser.add_argument('-out_molecule', type=str)
    
    args = parser.parse_args()
    
    print("Printing Mol")

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.out_molecule.split(".")[1])

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.molecule) 
    
    printResidueInfo(mol)
    
    phi_dihedral = np.zeros(mol.NumResidues(), dtype=np.float32)
    psi_dihedral = np.zeros(mol.NumResidues(), dtype=np.float32)
    
    for i in range (0, mol.NumResidues()):
        obres = mol.GetResidue(i)
        
        phi_dihedral = np.random.uniform(low=-180, high=180)
        psi_dihedral = np.random.uniform(low=-180, high=180)
        
        print(phi_dihedral, psi_dihedral)
        mi.SetPhiPsi(mol, obres, phi_dihedral, psi_dihedral)
    
    print("Writing file")
    obConversion.WriteFile(mol, args.out_molecule)
    
    
    
    
    
