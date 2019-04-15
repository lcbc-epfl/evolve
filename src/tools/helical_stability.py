'''
main.py

@author: Nicholas Browning
'''

import sys
import argparse 
import openbabel
from src import outprocesses as op
import subprocess
from src import MoleculeInfo as mi
from src import MoleculeCreator as mc

def amber_energy_minimize(mol):
     
    directory = '.'

    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(mol, directory + "/mol.pdb")

    op.runtleap(work_dir=directory + "/", mol_file='mol.pdb', tleaptemp='tleap_template.in', tleapin="leap.in")

    op.runPMEMD(work_dir=directory + "/", np=16, amberin='amber_minimization.in')
    
    fout = open(directory + '/min_struct.pdb', 'w')
    ferr = open(directory + '/min_struct.log', 'w')
    try:
        proc = subprocess.Popen(["ambpdb", "-p", directory + "/mol.prmtop", "-c", directory + "/amber.rst"], stdout=fout, stderr=ferr)
        proc.wait()
        fout.close()
        ferr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("convert rst to pdb failed, returned code %d (check '" + directory + "/min_struct.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to convert rst to pdb: %s" % (str(e)))
    
    f = open(directory + "/amber.out", "r")
    
    obConversion = openbabel.OBConversion() 
    obConversion.SetInFormat("pdb")

    molec = openbabel.OBMol()
    obConversion.ReadFile(molec, directory + '/min_struct.pdb')
  
    finalEnergy = op.parseAmberEnergy(directory + "/amber.out")
    
    return molec, finalEnergy

        
def helical_stability(initial_molecule, mutant_molecule, initialE, mutantE, mutation_points, dielectric):
    from src import constants
    from src import MoleculeInfo as mi
    directory = './'
               
    add = 0.0
    negate = initialE
    
    print "calcing: ", i, [mi.getResType(mutant_molecule.GetResidue(j)) for j in mutation_points]
    
    for j in xrange (0, len(mutation_points)):
        res = mi.getResType(mutant_molecule.mol.GetResidue(mutation_points[j]))

        add += constants.energies['ALA'][dielectric]  # TODO this won't work for a non poly ala protein!!
        negate += constants.energies[res][dielectric]
     
        print "min_energy:", mutantE, add, negate
        mutantE += (add - negate)
        
    return mutantE
        
    
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
    parser.add_argument('-im', '--initial_molecule', type=str)
    parser.add_argument('-m', '--molecule', type=str)
    parser.add_argument('-mi', '--mutate_indexes', nargs="*")
    parser.add_argument('-d', '--dielectric', type=int)
    args = parser.parse_args()
    
    print args

    obConversion = openbabel.OBConversion()
    
    obConversion.SetInAndOutFormats(args.initial_molecule.split(".")[1], args.initial_molecule.split(".")[1])

    initial_mol = openbabel.OBMol()
    

    if (not obConversion.ReadFile(initial_mol, args.initial_molecule)) :
        print "Couldn't read:", args.molecule
        sys.exit()
    print "Succesfully Read Molecule:", args.initial_molecule
    initial_mol_opt, initial_molE = amber_energy_minimize(initial_mol)
  
    print initial_molE
    
    mol = openbabel.OBMol()
    
    if (not obConversion.ReadFile(mol, args.molecule)) :
        print "Couldn't read:", args.molecule
        sys.exit()
    print "Succesfully Read Molecule:", args.molecule
    
    #mutant_mol_opt, mutant_molE = amber_energy_minimize(mol)
    
    stab = 0
    #stab = helical_stability(args.initial_mol_opt, args.mutant_mol_opt, initial_molE, mutant_molE, args.mutate_indexes, args.dielectric)

    print stab
    
