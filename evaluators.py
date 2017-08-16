'''
main.py

@author: Nicholas Browning
'''

from src.gaapi import Individual
import numpy as np
import openbabel


def testDihedralFitness(settings, individuals, fitness_index):
    # test with optimal being 180, -180 for phi, psi
    
    for i in xrange (0, settings.population_size):
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

    for i in xrange (0, settings.population_size):
        
        if (settings.composition_optimization):
            individuals[i].applyComposition(settings)

        if (settings.dihedral_optimization):
            individuals[i].applyPhiPsiDihedrals(settings)
        
        
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
         #
        if ff.Setup(individuals[i].mol) == 0:
            print "Could not setup forcefield"


        ff.SteepestDescent(500)

        ff.GetCoordinates(individuals[i].mol)
        
        phi_psi_dihedrals = individuals[i].getPhiPsiDihedrals()
            
        for j in xrange (0, len(phi_psi_dihedrals)):
            individuals[i].phi_dihedrals[j] = phi_psi_dihedrals[j][0]
            individuals[i].psi_dihedrals[j] = phi_psi_dihedrals[j][1]
            
        individuals[i].setFitness(fitness_index, ff.Energy(False)) 
        
        print i, individuals[i].fitnesses
        
def AmberBackboneStability(settings, individuals):
    import src.RotamerCreator as rc
    

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
        
        constraints.AddAtomConstraint(atom_idx)  # #add backbone atoms as constraints
            
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
    
    for i in xrange (0, settings.population_size):
        
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
        
        # $energy      SCF               SCFKIN            SCFPOT
        #     1 -350.2762671252      348.5663740952     -698.8426412204
        # $end

