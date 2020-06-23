'''
Evaluators in EVOLVE are functions that can call outprocesses and that return fitness values.



.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu

'''
from __future__ import print_function

from builtins import str
from builtins import range
from src.gaapi import Individual
from src.gaapi import generators
from src import constants as cnts
import numpy as np
import openbabel
import src.outprocesses as op
import os
import subprocess
import sys
import time
import glob
from src import MoleculeInfo as mi

def minimise_sidechain_ff(settings, individual):
    '''
    Use openbabel to minimize the sidechain

    Parameters
    ----------
    settings : object
            see :class:`src.JobConfig.Settings`
    individual: object
            see :class:`src.gaapi.Inidividual`


    '''

    from src import MoleculeInfo as mi
    
    constraints = openbabel.OBFFConstraints()
    
    backbone_atoms = []
    
    for i in range (0, individual.mol.NumResidues()):
        
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
    
    print("num_backbone_atoms", len(backbone_atoms))
    
    for atom_idx in backbone_atoms:
        
        constraints.AddAtomConstraint(atom_idx)  # add backbone atoms as constraints
            
    ff = openbabel.OBForceField.FindForceField("MMFF94")
        
    if ff == 0:
        print("Could not find forcefield")

    ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 

    ff.SetLogToStdErr() 
        
    if ff.Setup(individual.mol, constraints) == 0:
        print("Could not setup forcefield")

    ff.SteepestDescent(300)

    ff.GetCoordinates(individual.mol)


def getSidechainAngleAtoms(mol, sidechain_atoms):
    '''



    Parameters
    ----------
    mol
    sidechain_atoms

    Returns
    -------

    '''

    angle_atom_dx = []
    for obangle in openbabel.OBMolAngleIter(mol):
        
        vertexidx, atom1idx, atom2idx = obangle
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
    
    for i in range (0, individual.mol.NumResidues()):
        
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
        
    print(individual.mol.NumAtoms(), individual.mol.NumBonds())
    print("data:", len(bond_constraints), len(angle_constraints), len(dihedral_constraints))
    
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
        print("Could not find forcefield")

    ff.SetLogLevel(openbabel.OBFF_LOGLVL_NONE) 

    ff.SetLogToStdErr() 
        
    if ff.Setup(individual.mol, constraints) == 0:
        print("Could not setup forcefield")

    ff.SteepestDescent(300)

    ff.GetCoordinates(individual.mol)

    
def turbomol_scf_energy(settings, individuals, fitness_index):
    '''
    Does a turbomole SCF calculation

    Parameters
    ----------
    settings
    individuals
    fitness_index

    Returns
    -------

    '''
    import subprocess
    
    for i in range (0, len(individuals)):
        
        ind = individuals[i]
        
        minimise_backbone_ff(settings, ind)
        
        obConversion = openbabel.OBConversion()
        
        obConversion.SetOutFormat("xyz")
        
        obConversion.WriteFile(ind.mol, "coord.xyz")
        
        p = subprocess.Popen('x2t coord.xyz > coord', shell=True)
        
        p = subprocess.Popen('define < ' + settings.turbomole_template, shell=True)
        
        print("Starting SCF")
        
        with open('dscf.out', "w") as outfile:
            subprocess.call('dscf', stdout=outfile)
            
        print("Finished SCF")
        
        enfile = open('energy', 'r')
        
        lines = enfile.readlines()
        
        enfile.close()
        
        data = lines[1].split()
        
        ind.setFitness(fitness_index, float(data[1]))
        
        print(i, individuals[i].fitnesses)


def amber_energy_minimize(settings, individual):
    '''

    This function performs necessary topology preparation and minimization of a pdb file.

    The following evaluators depend on this:
    used in :meth:`helical_stability`

    Procedure:

    - write out mutated file from memory to pdb file on disk for amber
    - :meth:`src.outprocesses.runtleap`
    - :meth:`src.outprocesses.runPMEMD`
    - convert rst file to pdb file using `ambpdb`
    - read in the minimized structure as new pdb file and set as structure for  :data:`individual`

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individual : object
        see :class:`src.gaapi.Individual.Individual`

    Returns
    -------
    finalEnergy: float
        amber energy after minimization
    '''
    import shutil
     
    directory = settings.output_path + "/amber_run"
    
    if (os.path.exists(directory)):
        shutil.rmtree(directory, ignore_errors=True) 
        
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    if (os.path.exists(directory + "/mol.pdb")):
        os.remove(directory + "/mol.pdb")

    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(individual.mol, directory + "/mol.pdb")

    op.runtleap(work_dir=directory + "/", mol_file='mol.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in")

    op.runPMEMD(work_dir=directory + "/", np=settings.mpi_procs, amberin=settings.amber_params)
    
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
    individual.mol = molec
    
    finalEnergy = op.parseAmberEnergy(directory + "/amber.out")
    f.close()
    ferr.close()
    fout.close()
    return finalEnergy

def openmm_energy_minimize(settings, individual):
    '''

    This function performs necessary topology preparation and minimization of a pdb file.

    The following evaluators depend on this:
    used in :meth:`helical_stability`

    This uses OpenMM instead of Amber for the minimization. 
    There are options to either relax the structure using
    
    - a Steepest Descent Integrator
    - Simulated Annealing 
    - using minimization, MD (100ps) and then minimization

    Procedure:

    - write out mutated file from memory to pdb file on disk for amber
    - :meth:`src.outprocesses.runtleap`
    - openmm
    - read in the minimized structure as new pdb file and set as structure for  :data:`individual`

    .. code-block:: none

        [EVALUATOR]
        evaluators = helical_stability
        tleap_template = /path/to/tleap_template.in
        amber_params =  /path/to/amber_minimization.in
        dielectric = 80.0
        gpu_openmm = True
        gpudeviceindex = 0
        simAnneal = True 
        md = False

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individual : object
        see :class:`src.gaapi.Individual.Individual`

    Returns
    -------
    finalEnergy: float
        amber energy after minimization
    '''
    import shutil
    # OpenMM Imports
    import simtk.openmm as mm
    import simtk.openmm.app as app
    from parmed import load_file, unit as u
    from parmed.openmm import StateDataReporter, NetCDFReporter
    import openmmtools
    

    directory = settings.output_path + "/amber_run"
    
    if (os.path.exists(directory)):
        shutil.rmtree(directory, ignore_errors=True) 
        
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # if (os.path.exists(directory + "/mol.pdb")):
    #     os.remove(directory + "/mol.pdb")

    obconv = openbabel.OBConversion()
    obconv.SetOutFormat("pdb")
    obconv.WriteFile(individual.mol, directory + "/mol.pdb")

    op.runtleap(work_dir=directory + "/", mol_file='mol.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in")
    
    tleapfiles = load_file(os.getcwd()+"/amber_run/mol.prmtop",os.getcwd()+"/amber_run/mol.inpcrd")
    system = tleapfiles.createSystem(nonbondedMethod=app.NoCutoff,implicitSolvent=app.OBC2,)
   
    platform = mm.Platform.getPlatformByName('CUDA')
    prop = dict(CudaPrecision='mixed') # Use mixed single/double precision
    platform.setPropertyDefaultValue("CudaDeviceIndex", str(settings.gpudeviceindex))

    if settings.md:
        try:
            # Energy Minimize, then 100ps of NVT, then Gradient Descent Minimization.
            print("Using MD relaxation")
            integrator = mm.LangevinIntegrator(300*u.kelvin,1.0/u.picoseconds,2.0*u.femtoseconds,)
            sim = app.Simulation(tleapfiles.topology, system, integrator, platform, prop)
            sim.context.setPositions(tleapfiles.positions)
            sim.minimizeEnergy(maxIterations=2500)
            sim.reporters.append(StateDataReporter(sys.stdout, 100000, step=True, potentialEnergy=True,kineticEnergy=True, temperature=True))
            sim.step(50000)
            state = sim.context.getState(getEnergy=True,getPositions=True)
            # Gradient Descent to get to potential energy well. 
            integrator = openmmtools.integrators.GradientDescentMinimizationIntegrator(initial_step_size=0.04 * u.angstroms)  
            sim = app.Simulation(tleapfiles.topology, system, integrator, platform, prop)
            sim.context.setPositions(state.getPositions())
            sim.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True,kineticEnergy=True, temperature=True))
            sim.reporters.append(app.PDBReporter(directory +'/min.pdb', 10000))
            sim.step(15000)
            state = sim.context.getState(getEnergy=True)
            finalEnergy=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
            pass
        except ValueError:
            finalEnergy=999999
            pass
    else:
        integrator = openmmtools.integrators.GradientDescentMinimizationIntegrator(initial_step_size=0.04 * u.angstroms)  
        sim = app.Simulation(tleapfiles.topology, system, integrator, platform, prop)
        sim.context.setPositions(tleapfiles.positions)
        
        sim.reporters.append(StateDataReporter(os.devnull, 15000, step=True, potentialEnergy=True,kineticEnergy=True, temperature=True))
        sim.reporters.append(app.PDBReporter(directory +'/min.pdb', 15000))
        sim.step(15000)
        state = sim.context.getState(getEnergy=True,getPositions=True)
        finalEnergy=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
        
        obConversion = openbabel.OBConversion() 
        obConversion.SetInFormat("pdb")

        molec = openbabel.OBMol()
        obConversion.ReadFile(molec, directory + '/min.pdb')
        individual.mol = molec

    # ToDo:  
    # Implement settings for starting temperature and MD 
    if settings.simAnneal:
        try:
            minimized=state.getPositions()
            integrator = mm.LangevinIntegrator(300.0*u.kelvin,1.0/u.picoseconds,1.0*u.femtoseconds)
            sim = app.Simulation(tleapfiles.topology, system, integrator, platform, prop)
            sim.context.setPositions(minimized)
            sim.reporters.append(StateDataReporter(os.devnull, 2500, step=True, potentialEnergy=True,kineticEnergy=True, temperature=True))
            #sim.reporters.append(app.DCDReporter(name+'.dcd', 10000))
            for i in range(121):
                integrator.setTemperature(3*(120-i)*u.kelvin)
                sim.step(5000)
            state = sim.context.getState(getEnergy=True)
            finalEnergy=state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
            pass
        except ValueError:
            finalEnergy=999999
            pass
        

    return finalEnergy


def amber_energy_simplified(settings, individuals, fitness_index, pop_start=0):
    '''

    Remenant of Nick. Should probably be removed.

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individuals : object
        see :class:`src.gaapi.Individual.Individual`
    fitness_index
    pop_start

    Returns
    -------

    '''
     
    for i in range(pop_start, len(individuals)):
        
        print("Minimising: ", i, end=' ')
        
        finalEnergy = amber_energy_minimize(settings, individuals[i])
        print("min_energy:", finalEnergy)
        
        individuals[i].setFitness(fitness_index, finalEnergy)
        
        if (settings.solution_min_fitness is not None):
            if (individuals[i].getFitness(fitness_index) < settings.solution_min_fitness):
                individuals[i].setFitness(fitness_index, 999999999.0)
    
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))

        
def helical_stability(settings, individuals, fitness_index, pop_start=0):
    '''
    Helical stability function that computes the stability using an implicit solvent minimized structure.

    First the evaluator checks whether the current individual is identical in composition to an already computed individual, if yes computation is skipped and fitness value is copied.
    If no,  :meth:`amber_energy_minimize`

    If GPU option is set, openmm_minimize is used. 

    Then the add and negate energies and the inital energy will be added to the minimized energy

    To use this evaluator in the input file you need to include the following section in your file.

    .. code-block:: none

        [EVALUATOR]
        evaluators = helical_stability
        tleap_template = /path/to/tleap_template.in
        amber_params =  /path/to/amber_minimization.in
        mpi_processors = 20
        dielectric = 80.0
        gpu_openmm = False


    Reference: Perez et al., "EVOLVE: A Genetic Algorithm to Predict Protein Thermostability."

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individuals : list
        every item is an :class:`src.gaapi.Individual.Individual`
    fitness_index : int
        the index in the fitness list that helical_stability populates
    pop_start : int
        Defines with which individual to start


    Returns
    -------

    '''
    from src import constants
    from src import MoleculeInfo as mi
    directory = settings.output_path + "/amber_run"
     
    for i in range(pop_start, len(individuals)):
        
        already_done = -1
        
        if (i > 0):
            for j in range(0, i):
                if (np.array_equal(individuals[j].composition, individuals[i].composition)):
                    already_done = j
                    break
                    
        if (already_done != -1):
            print("Already computed: " , i, " -> member ", already_done)
            individuals[i].mol = openbabel.OBMol(individuals[already_done].mol)
            individuals[i].fitnesses = individuals[already_done].fitnesses
            continue     
               
        add = 0.0
        negate = settings.initial_energy
    
        print("Minimising: ", i, [mi.getResType(individuals[i].mol.GetResidue(j)) for j in settings.composition_residue_indexes])
        print("Rotamers: ", [cnts.selected_rotamers[v] for v in individuals[i].composition])
        
        if settings.gpu_openmm:
            finalEnergy = openmm_energy_minimize(settings, individuals[i])
        else:
            finalEnergy = amber_energy_minimize(settings, individuals[i])
        
        for j in range (0, len(settings.composition_residue_indexes)):
            #res = mi.getResType(individuals[i].mol.GetResidue(settings.composition_residue_indexes[j]))
            #add += constants.energies['ALA'][settings.helical_dielectric]  # TODO this won't work for a non poly ala protein!!
            #negate += constants.energies[res][settings.helical_dielectric]
            res = mi.getResType(individuals[i].mol.GetResidue(settings.composition_residue_indexes[j]))
            
            add += constants.energies[str(settings.originalResidues[j])][settings.helical_dielectric]
            negate += constants.energies[res][settings.helical_dielectric]
            print('add',str(settings.originalResidues[j]), constants.energies[str(settings.originalResidues[j])][settings.helical_dielectric])
            print('negate',res,'-(', settings.initial_energy,'+', constants.energies[res][settings.helical_dielectric],')')
            
        print("min_energy:", finalEnergy, add, negate)
        finalEnergy += (add - negate)
        
        individuals[i].setFitness(fitness_index, finalEnergy)
        
        if (settings.solution_min_fitness is not None):
            if (individuals[i].getFitness(fitness_index) < settings.solution_min_fitness):
                individuals[i].setFitness(fitness_index, 999999999.0)
                
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))




def stability_multi(settings, individuals, fitness_index, pop_start=0):
    '''

    This evaluator is used in conjunction with the :func:`evaluators.mmpbsa_multi`.

    It uses the total energy extracted from the mmpbsa calculation which is stored as the mean of all frames in `work_dir/evaluator_stability_multi.log` and the energy corrections developed by Perez et.al.
    to calculate the stability of a mutant
     Reference: Perez et al., "EVOLVE: A Genetic Algorithm to Predict Protein Thermostability.

    .. code-block:: none

        [GA_MULTIMMPBSA]
        multi_individual = True
        no_frames =12
        molecule_dir=./pdbs_rename/
        use_compute_cluster=False
        compute_cluster_nodes=6 #(optional)
        compute_cluster_ntasks=10 #(optional)
        compute_cluster_queuename=nccr-short #(optional)
        energy_calculator=pb
        [EVALUATOR]
        evaluators = mmpbsa_multi stability_multi
        tleap_template =path/to/leap.in
        amber_params =/path/to/min.in
        mmpbsa_params=/path/to/mmpbsa.in
        mpi_processors =24
        dielectric = 80.0

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individuals: list
        list of  :class:`src.gaapi.Individual.Individual`
    fitness_index: int
        integer designating the fitness index of this fitness (needed to update the correct value in the fitness list
    pop_start: int
        if only treating part of the population

    Returns
    -------

    '''
    from src import constants
    from src import MoleculeInfo as mi
    directory = settings.output_path + "/work_dir"

    for i in range(pop_start, len(individuals)):

        already_done = -1
        if (i > 0):
            for j in range(0, i):
                if (np.array_equal(individuals[j].composition, individuals[i].composition)):
                    already_done = j
                    break
        if (already_done != -1):
            print("Already computed: ", i, " -> member ", already_done)
            individuals[i].fitnesses = individuals[already_done].fitnesses
            continue
        if (individuals[i].fitnesses[0] == 123456789):
            # TODO set to max fitness in this generation to prevent replication
            individuals[i].setFitness(0, 3)
            individuals[i].setFitness(1, 100)
            continue

        add = 0.0
        negate = settings.initial_energy

        print("Getting energy for ", i, [mi.getResType(individuals[i].mol[0].GetResidue(j)) for j in settings.composition_residue_indexes])

        logfileStability = open(directory + "/evaluator_stability_multi.log")
        # lines=logfileStability.readlines()
        if individuals[i].stability == "NA":
            print("Stability=NA, stability set to 100")
            finalEnergy = 100
            pass
        else:
            finalEnergy = individuals[i].stability
            for j in range(0, len(settings.composition_residue_indexes)):
                res = mi.getResType(individuals[i].mol[0].GetResidue(settings.composition_residue_indexes[j]))
                add += constants.energies[str(settings.originalResidues[j])][settings.helical_dielectric]
                negate += constants.energies[res][settings.helical_dielectric]
                # print "Residue:", res, "add", add, "negate", negate
            print
            add, negate, finalEnergy
            finalEnergy += (add - negate)
        logfileStability.close()

        # if (i+1==len(individuals)):
        # truncate logfileStability after evaluation has finished
        #    logfileStability = open(directory+"/evaluator_stability_multi.log", "w")
        #    logfileStability.close()

        individuals[i].setFitness(fitness_index, finalEnergy)
        print("Fitnes index:", fitness_index)
        print("Fitness value:", finalEnergy)
        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))


def mmpbsa_multi(settings, individuals, fitness_index, pop_start=0):
    '''

    This evaluator is used in conjunction with the stability multi evaluator.

    It sets up the system (containing multiple frames) using leap, runs energy minimization (explicit solvent),
    then runs an MMPBSA calculation and extracts total energy (used for :func:`evaluators.stability_multi`) and binding free energy.

    .. code-block:: none

        [GA_MULTIMMPBSA]
        multi_individual = True
        no_frames =12
        molecule_dir=./pdbs_rename/
        use_compute_cluster=False
        compute_cluster_nodes=6 #(optional)
        compute_cluster_ntasks=10 #(optional)
        compute_cluster_queuename=nccr-short #(optional)
        energy_calculator=pb
        [EVALUATOR]
        evaluators = mmpbsa_multi stability_multi
        tleap_template =/path/to/leap.in
        amber_params =/path/to/min.in
        mmpbsa_params=/path/to/mmpbsa.in
        mpi_processors =24
        dielectric = 80.0

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individuals: list
        list of  :class:`src.gaapi.Individual.Individual`
    fitness_index: int
        integer designating the fitness index of this fitness (needed to update the correct value in the fitness list
    pop_start: int
        if only treating part of the population
    '''
    # create working directory for all different frames, this is needed to avoid output garbling from leap, pdb4amber and amber
    work_dir = settings.output_path + "/work_dir"
    logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
    # iterate over all invdividuals in the population
    for i in range(pop_start, len(individuals)):
        logfile.write("\n")
        logfile.write("#" * 100)
        logfile.write("\n")
        individualStart = time.time()
        # For each individual we determine whether it has already been computed in this population. Each combination of amino acids is only computed once per generation
        already_done = -1
        if (i > 0):
            for j in range(0, i):
                if (np.array_equal(individuals[j].composition, individuals[i].composition)):
                    already_done = j
                    break
        if (already_done != -1):
            print("Already computed: ", i, " -> member ", already_done)
            logfile.write("Already computed: " + str(i) + " -> member " + str(already_done))
            logfile.write("\n")
            # Not sure why this was included for helical stability, we initialize a new instance of OBMol for each frame anyway so this is not needed for mmpbsa_multi
            # individuals[i].mol = [None]*settings.no_frames
            # for x in range(0, settings.no_frames):
            #     individuals[i].mol[x] = openbabel.OBMol(individuals[already_done].mol[x])
            individuals[i].fitnesses = individuals[already_done].fitnesses
            # need to write the second fitness also to the stability logfile
            # logfileStability = open(work_dir+"/evaluator_stability_multi.log", "a")
            # logfileStability.write(str(individuals[already_done].fitnesses[1])+"\n")
            # logfileStability.close()
            continue
            # check if work dir exists
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        # create subdirectories for all frames specified and if they exists delete all files in them
        for frame in range(0, settings.no_frames):
            if not os.path.exists(work_dir + "/frame" + str(frame)):
                os.makedirs(work_dir + "/frame" + str(frame))
            if (os.path.exists(work_dir + "/frame" + str(frame) + "/mutated.pdb")):
                for f in glob.glob(work_dir + "/frame" + str(frame) + "/*"):
                    os.remove(f)
            pass
        pass
        # for each frame mutate to individual composition
        # note that this required changing Individual.py
        # With evaluator mmpbsa multi set, individual will not hold a single mol but multiple instances of OBmol
        for x in range(0, settings.no_frames):
            obconv = openbabel.OBConversion()
            obconv.SetOutFormat("pdb")
            obconv.WriteFile(individuals[i].mol[x], work_dir + "/frame" + str(x) + "/mutated.pdb")
            pass
        pass

        print("Individual ind {}".format(i))
        logfile.write("Individual ind {}".format(i))
        logfile.write("\n")
        logfile.write("Individual composition:" + str([mi.getResType(individuals[i].mol[0].GetResidue(a)) for a in settings.composition_residue_indexes]))
        logfile.write("\n")
        logfile.close()
        # running tleap for all frames. Inside this routine we change to every frame dir and run a seperate instance of leap
        op.runtleapMultiInd(work_dir=work_dir, leap_inputfile=settings.tleap_template, n_frames=settings.no_frames)
        # run pmemd using multipmemd version. Need to use pmemd.MPI because of force truncation in pmemd.cuda!
        minReturnCode = op.runPMEMDMultiInd(work_dir=work_dir, np=settings.mpi_procs, amberin=settings.amber_params,
                                            n_frames=settings.no_frames, compute_cluster=settings.use_compute_cluster,
                                            nodes=settings.compute_cluster_nodes,
                                            ntasks=settings.compute_cluster_ntasks,
                                            queuename=settings.compute_cluster_queuename)
        logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
        if minReturnCode == 0:
            print("Setting indivdual fitness to 0 due to error in the minimization (Seg fault)")
            finalEnergy = 0
            pass
        else:
            logfile.write("finished energy min")
            logfile.write("\n")
            logfile.close()
            # Do mmpbsa calculation after frames have been minimized and converted to a solvent free pdb file. MMPBSA calculation can run on as many cores as there are frames.
            finalEnergy = op.runMMPBSAMultiInd(work_dir=work_dir, mmpbsa_inputfile=settings.mmpbsa_params,
                                               n_frames=settings.no_frames,
                                               compute_cluster=settings.use_compute_cluster,
                                               nodes=settings.compute_cluster_nodes,
                                               ntasks=settings.compute_cluster_ntasks,
                                               queuename=settings.compute_cluster_queuename,
                                               energy_evaluator=settings.energy_calculator)
            logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
            logfile.write("finished MMPBSA min\n")

            # set fitness to absurd value in case something goes wrong in order to avoid replication of this composition
            if (settings.solution_min_fitness is not None):
                if (finalEnergy < settings.solution_min_fitness):
                    finalEnergy = 0
        # update fitness of individual with computed value from MMPBSA GB or PB calculatioN
        individuals[i].setFitness(fitness_index, finalEnergy)
        if fitness_index == 0:
            individuals[i].setFitness(fitness_index + 1, 0)
            pass
        # read out energies from minimization and store in logfile in case evaluator stability multi is also used
        frameEnergies = settings.no_frames * [None]
        for frame in range(0, settings.no_frames):
            # frameEnergies[frame]=op.parseAmberEnergy(work_dir+"/frame"+str(frame)+"/min.out."+"frame"+str(frame))
            frameEnergies[frame] = op.parseGBEnergy(work_dir + "/_MMPBSA_complex_gb.mdout." + str(frame))
            logfile.write("frame" + str(frame) + " energy:" + str(frameEnergies[frame]) + "\n")
            if frameEnergies[frame] == 999:
                print("Something went wrong with reading the minimization file")
                continue
            pass
        if frameEnergies[frame] != 999:
            stabilityEnergy = np.mean(frameEnergies)
            logfile.write("Stability Energy (mean of all frames):" + str(stabilityEnergy) + "\n")
            logfile.write("Stability Energy  - Initial Energy:" + str(stabilityEnergy - settings.initial_energy) + "\n")
        else:
            stabilityEnergy = "NA"

        logfileStability = open(work_dir+"/evaluator_stability_multi.log", "a")
        logfileStability.write(str(stabilityEnergy)+"\n")
        logfileStability.close()
        individuals[i].setStability(stabilityEnergy)

        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))
        logfile.write('ind {}, fitness {}\n'.format(i, individuals[i].fitnesses))
        individualEnd = time.time()
        # read in new files
        # update mols with minized molecules from work dir

        if settings.update_files == True:
            for frame in range(0, settings.no_frames):
                errorcatcher = 0
                obconv = openbabel.OBConversion()
                pdbFile = openbabel.OBMol()
                obconv.ReadFile(pdbFile, settings.output_path + "/work_dir/frame" + str(frame) + "/min_nosolv.pdb")
                # TODO: this should only print once
                try:
                    # test if pdb file is correct
                    testPDBfile = [mi.getResType(pdbFile.GetResidue(j)) for j in settings.composition_residue_indexes]
                    pass
                except Exception as e:
                    print
                    "could not read pdb file", str(e)
                    individuals[i].setFitness(0, 123456789)
                    individuals[i].setFitness(1, 123456789)
                    errorcatcher = 1
                    pass
                # if minimization fails we cannot update with non exisiting pdb file
                if errorcatcher == 0:
                    individuals[i].mol[frame] = pdbFile
                    pass
                pass
        logfile.write("Runtime of this ind: " + str(individualEnd - individualStart) + " sec.")
        logfile.write("\n")
        logfile.write("#" * 100)
        logfile.write("\n")
        if minReturnCode != 0:
            os.remove(
                work_dir + "/reference.frc")  # This file is generated by MMPBSA and gets huge. Somehow MMPBSA does not remove it automatically
            # Get a list of all the files with _MMPBSA
            fileList = glob.glob(work_dir + "/_MMPBSA*")
            # Iterate over the list of filepaths & remove each file.
            for filePath in fileList:
                os.remove(filePath)
    pass
    logfile.close()
    
