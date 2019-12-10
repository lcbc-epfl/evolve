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


def printResidueInfo(mol):
    print("Num Atoms:", mol.NumAtoms())
    print("Num Bonds:", mol.NumBonds())
    print("Num Residues:", mol.NumResidues())
    
    
    print("Phi, Psi Dihedrals:")
    dihedrals = mi.getAllPhiPsiDihedrals(mol)

    (sorted_atoms_ID, sorted_atoms_OB) = mi.get_atoms_per_residue(mol)
    
    for i in range (0, mol.NumResidues()):
        res = mol.GetResidue(i)
        
        print("Res", i, res.GetName())

        print("Atoms {}".format(sorted_atoms_ID[i]))
        
        print("[phi, psi]", dihedrals[i][0], dihedrals[i][1])
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    
    args = parser.parse_args()
    
    print("Printing Mol")

    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat(args.molecule.split(".")[1])

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.molecule) 
    
    printResidueInfo(mol)
    
    
    
    
    
