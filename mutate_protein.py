'''
@author: Nicholas Browning
'''

# module load anaconda/4.3.1
# module load openbabel/2.3.2/gcc-4.9.2
# export PYTHONPATH=$PYTHONPATH:/path/to/GABio/src
# CALL ME LIKE THIS: python mutate_protein.py -molecule gpgg.pdb -rotamer_path /lcbcdata/nbrownin/Rotamers_mol2/LYS/K14.mol2 -out_molecule test.pdb -res_index 1


import argparse
import openbabel
from src import MoleculeCreator as mc
from src import MoleculeInfo as mi
import sys

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    parser.add_argument('-rotamer_path', type=str)
    parser.add_argument('-out_molecule', type=str)
    parser.add_argument('-res_index', type=int)
    
    args = parser.parse_args()
    
    print "Printing Mol"

    obConversion = openbabel.OBConversion()
    
    # "test.xyz" -> ['test', 'xyz']
    
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.out_molecule.split(".")[1])

    mol = openbabel.OBMol()
    
    if (not obConversion.ReadFile(mol, args.molecule)) :
        print "Couldn't read:", args.molecule
        sys.exit()
    print "Succesfully Read Molecule:", args.molecule
    
    print args.rotamer_path.split(".")[1]
    frag = openbabel.OBMol()
    obConversion.SetInFormat(args.rotamer_path.split(".")[1])
    if (not obConversion.ReadFile(frag, args.rotamer_path)):
        sys.exit()
        print "Couldn't read:", args.rotamer_path
    
    curr = mol.GetResidue(args.res_index)
    
    mol_CA = mi.getAlphaCarbon(curr)
    mol_CB = mi.getBetaAtom(curr)
    mol_N = mi.getBBNitrogen(curr)
    mol_CO = mi.getBBCarboxyl(curr)
    
    if (mol_CA is None or mol_CB is None or mol_N is None or mol_CO is None):
        print "Could not find backbone atoms"
        print "CA:", mol_CA, "CB:", mol_CB, "N:", mol_N, "C(=O):", mol_CO
        sys.exit()
    
    print "Swapping Side Chain"
    mc.swapsidechain(mol, args.res_index, frag)
    print "Writing File:", args.out_molecule
    obConversion.WriteFile(mol, args.out_molecule)
    

