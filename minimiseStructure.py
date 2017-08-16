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
        
    ff = openbabel.OBForceField.FindForceField("MMFF94")
        
    if ff == 0:
        print "Could not find forcefield"

    ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 

    ff.SetLogToStdErr() 

#
    if ff.Setup(mol) == 0:
        print "Could not setup forcefield"


    ff.SteepestDescent(50000)

    ff.GetCoordinates(mol)
        
    print ff.Energy(False) 
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    parser.add_argument('-out_molecule', type=str)
    
    args = parser.parse_args()
    
    print "Printing Mol"

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.out_molecule.split(".")[1])

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, args.molecule) 
    
    printResidueInfo(mol)
    
    obConversion.WriteFile(mol, args.out_molecule)
    
    
    
    
    
    
    
