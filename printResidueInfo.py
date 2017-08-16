'''
printMoleculeInfo

@author: Nicholas Browning
'''
# Use this to print residue information
import argparse
from src import MoleculeInfo as mi
import openbabel


def printResidueInfo(mol):
    print "Num Atoms:", mol.NumAtoms()
    print "Num Bonds:", mol.NumBonds()
    print "Num Residues:", mol.NumResidues()
    
    
    print "Phi, Psi Dihedrals:"
    dihedrals = mi.getAllPhiPsiDihedrals(mol)
    
    for i in xrange (0, mol.NumResidues()):
        res = mol.GetResidue(i)
        
        print "Res", i, res.GetName()
        
        print "[phi, psi]", dihedrals[i][0], dihedrals[i][1]
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    
    args = parser.parse_args()
    
    print "Printing Mol"

    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat(args.molecule.split(".")[1])

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.molecule) 
    
    printResidueInfo(mol)
    
    
    
    
    
