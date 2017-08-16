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
    
    pointers = []
    pointer_names = []
    indexes = []
    for i in xrange (0, mol.NumResidues()):
        res = mol.GetResidue(i)
        
        res_name = mi.getResType(res)
        
        if (res_name == "PRO"):
            pointer_names.append("PRO")
            pointers.append(3)
        
        
        found_prepro = False
        
        if (i + 1 < mol.NumResidues()):
            
            res_1 = mol.GetResidue(i + 1)
            res_1_name = mi.getResType(res_1)
            
            if (res_1_name == "PRO"):
                 pointer_names.append("PREPRO")
                 pointers.append(2)
                 found_prepro = True
                 
        if (not found_prepro):
            if (res_name == "GLY"):
                pointer_names.append("GLY")
                pointers.append(1)
            else:
                pointer_names.append("GEN")
                pointers.append(0)
                
        indexes.append(i)
    print " ".join([str(v) for v in indexes])
    print " ".join([str(v) for v in pointers])
    print pointer_names
                
    
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
    
    
    
    
    
