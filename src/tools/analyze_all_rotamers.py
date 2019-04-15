'''
main.py

@author: Nicholas Browning
'''

import os
import sys
import argparse 
import openbabel
from src import outprocesses as op
import subprocess
from src import MoleculeInfo as mi
from src import MoleculeCreator as mc
from src import constants as cnts    

def mutate(mol, mutate_indexes, mutate_residues, rotamer_path='/share/lcbcsrv5/lcbcdata/nbrownin/Rotamer_Library_Flat'):
    
    for i in xrange(0, len(mutate_indexes)):

        frag = openbabel.OBMol()
        obConversion.SetInFormat(mutate_residues[i].split(".")[1])
        if (not obConversion.ReadFile(frag, rotamer_path + "/" + mutate_residues[i])):
            sys.exit()
            print "Couldn't read:", args.rotamer_path
    
        curr = mol.GetResidue(int(mutate_indexes[i]))
    
        mol_CA = mi.getAlphaCarbon(curr)
        mol_CB = mi.getBetaAtom(curr)
        mol_N = mi.getBBNitrogen(curr)
        mol_CO = mi.getBBCarboxyl(curr)
    
        if (mol_CA is None or mol_CB is None or mol_N is None or mol_CO is None):
            print "Could not find backbone atoms"
            print "CA:", mol_CA, "CB:", mol_CB, "N:", mol_N, "C(=O):", mol_CO
            sys.exit()
            
        mc.swapsidechain(mol, int(mutate_indexes[i]), frag)  
        
def sander_minimize(mol, num_procs):
    delete_files = ["mol.pdb", "mol.prmtop", "mol.inpcrd", "mol.rst", "leap.log", "amber.out", "amber.log", "err", "mdinfo", "mdout"]
    for v in delete_files:
        if os.path.exists(v):
            os.remove(v)
            
    obconv = openbabel.OBConversion()
    
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(mol, "mol.pdb")

    op.runtleap(work_dir="", mol_file='mol.pdb', tleaptemp="tleap_template.in", tleapin="leap.in")
    
    op.runAmberMPI(np=num_procs)
   
    finalEnergy = op.parseAmberEnergy("amber.out")
    
    return finalEnergy


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--molecule', type=str)
    parser.add_argument('-mi', '--mutate_indexes', nargs="*")
    parser.add_argument('-np', '--num_procs')
    
    args = parser.parse_args()

    obConversion = openbabel.OBConversion()
    
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.molecule.split(".")[1])
    
    mol = openbabel.OBMol()

    rotamers = cnts.rotamers
    
    if (not obConversion.ReadFile(mol, args.molecule)) :
        print "Couldn't read:", args.molecule
        sys.exit()
        
    print rotamers
    for i, rotamer in enumerate(rotamers):
        initial_mol = openbabel.OBMol(mol)
        print 'Calculating Energy: ', rotamer
        rotamer_list = [None] * len(args.mutate_indexes)
        for i in xrange(0, len(args.mutate_indexes)):
            rotamer_list[i] = rotamer + '.mol2'
        
        mutate(mol, args.mutate_indexes, rotamer_list, rotamer_path='/share/lcbcsrv5/lcbcdata/nbrownin/Rotamer_Library_Flat')
    
        initial_e = sander_minimize(mol, args.num_procs)
    
        print initial_e
    
         
        

    
    
