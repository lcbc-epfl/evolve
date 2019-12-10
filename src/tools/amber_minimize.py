'''
@author: Nicholas Browning
'''
from __future__ import print_function
from builtins import str
import os
import sys
import argparse 
import openbabel
from src import outprocesses as op
import subprocess


def sander_minimize(mol, out_mol, num_procs):
    delete_files = ["mol.pdb", "mol.prmtop", "mol.inpcrd", "mol.rst", "leap.log", "amber.out", "amber.log", "err", "mdinfo", "mdout"]
    for v in delete_files:
        if os.path.exists(v):
            os.remove(v)
            
    obconv = openbabel.OBConversion()
    
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(mol, "mol.pdb")

    op.runtleap(work_dir="", mol_file='mol.pdb', tleaptemp="tleap_template.in", tleapin="leap.in")
    
    op.runAmberMPI(np=num_procs)

    fout = open(out_mol, 'w')
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
    
    return finalEnergy


def pmemd_minimize(mol, out_mol, cuda, num_procs):

    delete_files = ["mol.pdb", "mol.prmtop", "mol.inpcrd", "mol.rst", "leap.log", "amber.out", "amber.log", "err", "mdinfo", "mdout"]
    for v in delete_files:
        if os.path.exists(v):
            os.remove(v)
    
    obconv = openbabel.OBConversion()
  
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(mol, "mol.pdb")

    op.runtleap(work_dir="", mol_file='mol.pdb', tleaptemp="tleap_template.in", tleapin="leap.in")
     
    op.runPMEMD(np=num_procs, use_cuda=cuda)
  
    fout = open(out_mol, 'w')
    ferr = open('err', 'w')
    try:
        proc = subprocess.Popen(["ambpdb", "-p", "mol.prmtop", "-c", "mol.rst"], stdout=fout, stderr=ferr)
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
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(out_mol.split(".")[1], out_mol.split(".")[1])

    obConversion.ReadFile(mol, args.molecule) 
    return finalEnergy
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-molecule', type=str)
    parser.add_argument('-cuda', type=bool)
    parser.add_argument('-out_molecule', type=str)
    parser.add_argument('-num_procs', type=str)
    parser.add_argument('-code', type=str)
    args = parser.parse_args()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(args.molecule.split(".")[1], args.out_molecule.split(".")[1])

    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, args.molecule) 
    
    if (args.code == "pmemd"):
        e = pmemd_minimize(mol, args.out_molecule, args.cuda, args.num_procs)
    else:
        e = sander_minimize(mol, args.out_molecule, args.num_procs)
    
    obConversion.WriteFile(mol, args.out_molecule) 
    print(e)
    
