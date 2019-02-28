'''
main.py

@author: Nicholas Browning
'''

from src.gaapi import Individual
from src.gaapi import generators
import numpy as np
import openbabel
import src.outprocesses as op
import os
import subprocess
import sys


def testDihedralFitness(settings, individuals, fitness_index):
    # test with optimal being 180, -180 for phi, psi
    
    for i in xrange (0, len(individuals)):
        individuals[i].applyPhiPsiDihedrals()
        fitness = 0.0
        phi, psi = individuals[i].phi_dihedrals, individuals[i].psi_dihedrals
        
        for j in xrange (0, len(settings.dihedral_residue_indexes)):
            
            if (phi[j] != 999):
                fitness += np.sqrt((180 - phi[j]) ** 2)
            if (psi[j] != 999):
                fitness += np.sqrt((180 - psi[j]) ** 2)
                
        individuals[i].setFitness(fitness_index, fitness)


def testEnergyByMinimisation(settings, individuals, fitness_index):

    # if (settings.verbose):
        # print("IND NUM, IND FITNESS")

    for i in xrange(0, len(individuals)):
        
        if (settings.composition_optimization and not individuals[i].init):
            individuals[i].applyComposition(settings)

        if (settings.dihedral_optimization and not individuals[i].init):
            individuals[i].applyPhiPsiDihedrals()

        if (settings.basilisk_and_sidechains and not individuals[i].init):
            individuals[i].apply_chi_dihedrals()

            
        # Minimize the structure in individuals[i].mol
        ff = openbabel.OBForceField.FindForceField("MMFF94")
        
        if ff == 0:
            print "Could not find forcefield"

        ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 
        ff.SetLogToStdErr()
        
#        # // Set the constraints
      # OBFFConstraints constraints;
      # constraints.AddAtomConstraint(1);

      # // We pass the constraints as argument for Setup()
      # if (!pFF->Setup(mol, constraints)) {
      # cerr << "ERROR: could not setup force field." << endl;
      # }
         #
        if ff.Setup(individuals[i].mol) == 0:
            print "Could not setup forcefield"


        ff.SteepestDescent(500)  # Perform the actual minimization, maximum 500 steps
        # ff.ConjugateGradients(500) #If using minimizing using conjugate gradient method (symmetric, positive definite matrix or quadratic functions)

        ff.GetCoordinates(individuals[i].mol)
        
        phi_psi_dihedrals = individuals[i].getPhiPsiDihedrals()
            
        for j in xrange (0, len(phi_psi_dihedrals)):
            individuals[i].phi_dihedrals[j] = phi_psi_dihedrals[j][0]
            individuals[i].psi_dihedrals[j] = phi_psi_dihedrals[j][1]

        if (settings.basilisk_and_sidechains):
            chi_dihedrals = individuals[i].get_chi_dihedrals_per_res()
            individuals[i].chi_angles = chi_dihedrals
            
        individuals[i].setFitness(fitness_index, ff.Energy(False))

        # if (settings.verbose):2>&1 | tee -i
           # print i, individuals[i].fitnesses
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))

def minimise_sidechain_ff(settings, individual):
    from src import MoleculeInfo as mi
    
    constraints = openbabel.OBFFConstraints()
    
    backbone_atoms = []
    
    for i in xrange (0, individual.mol.NumResidues()):
        
        res = individual.mol.GetResidue(i)
        alpha_carbon = mi.getAlphaCarbon(res)
        beta_atom = mi.getBetaAtom(res)
        
        if (mi.getResType(res) == "PRO"):
            # ignore PRO for now
            continue
        
        sidechain_atoms = openbabel.vectorInt()
        
        individual.mol.FindChildren(sidechain_atoms, alpha_carbon.GetIdx(), beta_atom.GetIdx())
        
        for curr_atom in openbabel.OBResidueAtomIter(res):
            
            if (curr_atom.GetIdx() in sidechain_atoms):
                continue
            
            backbone_atoms.append(curr_atom.GetIdx())
    
    print "num_backbone_atoms", len(backbone_atoms)
    
    for atom_idx in backbone_atoms:
        
        constraints.AddAtomConstraint(atom_idx)  # add backbone atoms as constraints
            
    ff = openbabel.OBForceField.FindForceField("MMFF94")
        
    if ff == 0:
        print "Could not find forcefield"

    ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 

    ff.SetLogToStdErr() 
        
    if ff.Setup(individual.mol, constraints) == 0:
        print "Could not setup forcefield"

    ff.SteepestDescent(300)

    ff.GetCoordinates(individual.mol)


def getSidechainAngleAtoms(mol, sidechain_atoms):

    angle_atom_dx = []
    for obangle in openbabel.OBMolAngleIter(mol):
        
        vertexidx, atom1idx, atom2idx = obangle
        
        # vertexidx = vertex.GetIdx()
        # atom1idx = atom1.GetIdx()
        # atom2idx = atom2.GetIdx()
        
        if (vertexidx in sidechain_atoms or atom1idx in sidechain_atoms or atom2idx in sidechain_atoms):
            angle = mol.GetAngle(mol.GetAtom(atom1idx), mol.GetAtom(vertexidx), mol.GetAtom(atom2idx))
            angle_atom_dx.append([atom1idx, vertexidx, atom2idx, angle])
    return angle_atom_dx
        

def getSidechainBondAtoms(mol, sidechain_atoms):
    bond_atom_dx = []
    for obbond in openbabel.OBMolBondIter(mol):
        
        atom1 = obbond.GetBeginAtom()
        atom2 = obbond.GetEndAtom()
    
        atom1idx = atom1.GetIdx()
        atom2idx = atom2.GetIdx()
        
        if (atom1idx in sidechain_atoms or atom2idx in sidechain_atoms):
            bondL = obbond.GetLength()
            bond_atom_dx.append([atom1idx, atom2idx, bondL])
            
    return bond_atom_dx


def getSidechainDihedralAtoms(mol, sidechain_atoms):
    torsion_atom_idx = []
    for obtorsion in openbabel.OBMolTorsionIter(mol):
        
        atom1idx, atom2idx, atom3idx, atom4idx = obtorsion
        

        if (atom1idx in sidechain_atoms or atom2idx in sidechain_atoms or atom3idx in sidechain_atoms or atom4idx in sidechain_atoms):
            torsion = mol.GetTorsion(atom1idx, atom2idx, atom3idx, atom4idx)
            torsion_atom_idx.append([atom1idx, atom2idx, atom3idx, atom4idx, torsion])

    return torsion_atom_idx


def minimise_backbone_ff(settings, individual):
    from src import MoleculeInfo as mi
    
    constraints = openbabel.OBFFConstraints()
    
    angle_constraints = []
    bond_constraints = []
    dihedral_constraints = []
    
    for i in xrange (0, individual.mol.NumResidues()):
        
        res = individual.mol.GetResidue(i)
        alpha_carbon = mi.getAlphaCarbon(res)
        beta_atom = mi.getBetaAtom(res)

        if (mi.getResType(res) == "PRO"):
            # ignore PRO for now
            continue
        
        ss_atoms = openbabel.vectorInt()
        
        individual.mol.FindChildren(ss_atoms, alpha_carbon.GetIdx(), beta_atom.GetIdx())
        
        ss_atoms = [int(atom) for atom in ss_atoms]
        
        ss_atoms.append(beta_atom.GetIdx())
        
        angle_info = getSidechainAngleAtoms(individual.mol, ss_atoms)
        bond_info = getSidechainBondAtoms(individual.mol, ss_atoms)
        torsion_info = getSidechainDihedralAtoms(individual.mol, ss_atoms)
        
        for v in bond_info:
            bond_constraints.append(v)
        for v in angle_info:
            angle_constraints.append(v)
        for v in torsion_info:
            dihedral_constraints.append(v)
        
    print individual.mol.NumAtoms(), individual.mol.NumBonds()
    print "data:", len(bond_constraints), len(angle_constraints), len(dihedral_constraints)
    
    # OBMolBondIter, OBMolAngleIter, OBMolTorsionIter, OBMolRingIter -
    # void     AddDistanceConstraint (int a, int b, double length)
    # void     AddAngleConstraint (int a, int b, int c, double angle)
    # void     AddTorsionConstraint (int a, int b, int c, int d, double torsion)
    
    for bond_constraint in bond_constraints:
        constraints.AddDistanceConstraint(bond_constraint[0], bond_constraint[1], bond_constraint[2])
    
    for ang_constraint in angle_constraints:
        constraints.AddAngleConstraint(ang_constraint[0], ang_constraint[1], ang_constraint[2], ang_constraint[3]) 
        
    for tor_constraint in dihedral_constraints:
        constraints.AddTorsionConstraint(tor_constraint[0], tor_constraint[1], tor_constraint[2], tor_constraint[3], tor_constraint[4]) 
            
    ff = openbabel.OBForceField.FindForceField("MMFF94")
        
    if ff == 0:
        print "Could not find forcefield"

    ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 

    ff.SetLogToStdErr() 
        
    if ff.Setup(individual.mol, constraints) == 0:
        print "Could not setup forcefield"


    ff.SteepestDescent(300)

    ff.GetCoordinates(individual.mol)

    
def turbomol_scf_energy(settings, individuals, fitness_index): 
    import subprocess
    
    for i in xrange (0, len(individuals)):
        
        ind = individuals[i]
        
        minimise_backbone_ff(settings, ind)
        
        obConversion = openbabel.OBConversion()
        
        obConversion.SetOutFormat("xyz")
        
        obConversion.WriteFile(ind.mol, "coord.xyz")
        
        p = subprocess.Popen('x2t coord.xyz > coord', shell=True)
        
        p = subprocess.Popen('define < ' + settings.turbomole_template, shell=True)
        
        print "Starting SCF"
        
        with open('dscf.out', "w") as outfile:
            subprocess.call('dscf', stdout=outfile)
            
        print "Finished SCF"
        
        enfile = open('energy', 'r')
        
        lines = enfile.readlines()
        
        enfile.close()
        
        data = lines[1].split()
        
        ind.setFitness(fitness_index, float(data[1]))
        
        print i, individuals[i].fitnesses


def amber_energy_minimize(settings, individual):
     
    directory = settings.output_path + "/amber_run"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    if (os.path.exists(directory + "/mol.pdb")):
        os.remove(directory + "/mol.pdb")

    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(individual.mol, directory + "/mol.pdb")

    op.runtleap(work_dir=directory + "/", mol_file='mol.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in")

    op.runPMEMD(work_dir=directory + "/", np=settings.mpi_procs, amberin=settings.amber_params)
    
    print "finished energy min"
    
    fout = open(directory + '/min_struct.pdb', 'w')
    ferr = open(directory + '/min_struct.log', 'w')
    try:
        proc = subprocess.Popen(["ambpdb","-p", directory + "/mol.prmtop", "-c", directory + "/amber.rst"], stdout=fout, stderr=ferr)
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
    individual.mol = molec
    
    finalEnergy = op.parseAmberEnergy(directory + "/amber.out")
    
    return finalEnergy

def amber_energy_simplified(settings, individuals, fitness_index, pop_start=0):
    
     
    for i in xrange(pop_start, len(individuals)):
        
        print "Minimising: ", i,
        
        finalEnergy = amber_energy_minimize(settings, individuals[i])
        
        individuals[i].setFitness(fitness_index, finalEnergy)
        
        if (settings.solution_min_fitness is not None):
            if (individuals[i].getFitness(fitness_index) < setings.solution_min_fitness):
                individuals[i].setFitness(fitness_index, 999999999.0)
        
    
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))
        
def helical_stability(settings, individuals, fitness_index, pop_start=0):
    from src import constants
    from src import MoleculeInfo as mi
    directory = settings.output_path + "/amber_run"
     
    for i in xrange(pop_start, len(individuals)):
        
        already_done = -1
        
        if (i > 0):
            for j in xrange(0, i):
                if (np.array_equal(individuals[j].composition, individuals[i].composition)):
                    already_done = j
                    break
                    
        if (already_done != -1):
            print "ALready computed: " , i, " -> member ", already_done 
            individuals[i].mol = openbabel.OBMol(individuals[already_done].mol)
            individuals[i].fitnesses = individuals[already_done].fitnesses
            continue     
               
        add = 0.0
        negate = settings.initial_energy
    
        print "Minimising: ", i, [mi.getResType(individuals[i].mol.GetResidue(j)) for j in settings.composition_residue_indexes]
        
        finalEnergy = amber_energy_minimize(settings, individuals[i])
        
        for j in xrange (0, len(settings.composition_residue_indexes)):
            res = mi.getResType(individuals[i].mol.GetResidue(settings.composition_residue_indexes[j]))
            add += constants.energies['ALA'][settings.helical_dielectric]
            negate += constants.energies[res][settings.helical_dielectric]
            
        print add, negate
        finalEnergy += (add - negate)
        
        individuals[i].setFitness(fitness_index, finalEnergy)
        
    
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))
        
    
    
    
    
