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

def pmemd_minimize(mol):

    obconv = openbabel.OBConversion()
    
    delete_files = ['mol.pdb', 'leap.in', 'amber.out', 'mol.prmtop', 'mol.inpcrd', 'mol.rst', 'mdinfo', 'amber.log', 'err', 'leap.log', 'logfile', 'amber.rst' ]
    for v in delete_files:
        if (os.path.exists(v)):
            os.remove(v)
    
    obconv.SetInAndOutFormats("pdb", "pdb")
    
    obconv.WriteFile(mol, "mol.pdb")

    op.runtleap(work_dir="", mol_file='mol.pdb', tleaptemp="tleap_template.in", tleapin="leap.in")
     
    op.runPMEMD(np=16)
  
    fout = open('tmp.pdb', 'w')
    ferr = open('err', 'w')
    try:
        proc = subprocess.Popen(["ambpdb", "-p", "mol.prmtop", "-c", "amber.rst"], stdout=fout, stderr=ferr)
        proc.wait()
        fout.close()
        ferr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("convert rst to pdb failed, returned code %d (check '" + "/min_struct.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to convert rst to pdb: %s" % (str(e)))

    
    finalEnergy = op.parseAmberEnergy("amber.out")    
    obconv.ReadFile(mol, 'tmp.pdb')
    
    return finalEnergy

def helical_stability(initial_mol, mol, mutated_indexes, DIE):
    from src import constants
    from src import MoleculeInfo as mi
     
    add = 0.0
    negate = 0.0
        
    for j in xrange (0, len(mutated_indexes)):
        res = mi.getResType(mol.GetResidue(int(mutated_indexes[j])))
        add += constants.energies['ALA'][DIE]
        negate += constants.energies[res][DIE]
   
    print add, negate
    return (add - negate)
        

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-m', '--molecule', type=str)
    parser.add_argument('-om', '--out_molecule', type=str)
    parser.add_argument('-mr', '--mutate_residues', nargs="*")
    parser.add_argument('-mi', '--mutate_indexes', nargs="*")
    parser.add_argument('-d', '--dielectric', type=int)
    args = parser.parse_args()
    

    obConversion = openbabel.OBConversion()
    
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.out_molecule.split(".")[1])

    mol = openbabel.OBMol()
    
    if (not obConversion.ReadFile(mol, args.molecule)) :
        print "Couldn't read:", args.molecule
        sys.exit()    
    
    initial_mol = openbabel.OBMol(mol)
    initial_e = pmemd_minimize(initial_mol)
    
    obConversion.WriteFile(initial_mol, 'initial_mol.pdb')
    
    mutate(mol, args.mutate_indexes, args.mutate_residues, rotamer_path='/share/lcbcsrv5/lcbcdata/nbrownin/Rotamer_Library_Flat')
    
    
    opt_e = pmemd_minimize(mol)
    
    diff = helical_stability(initial_mol, mol, args.mutate_indexes, args.dielectric)
    
    obConversion.WriteFile(mol, args.out_molecule)

    print opt_e - initial_e + diff, initial_e, opt_e, diff
