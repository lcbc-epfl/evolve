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


def amber_energy_simplified(settings, individuals, fitness_index, pop_start=0):
    for i in xrange(pop_start, len(individuals)):
        
        print "Minimising: ", i
  
        directory = settings.output_path + "/amber_run"
        
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        if (os.path.exists(directory + "/mol.pdb")):
            os.remove(directory + "/mol.pdb")

        obconv = openbabel.OBConversion()
        obconv.SetOutFormat("pdb")
        obconv.WriteFile(individuals[i].mol, directory + "/mol.pdb")
    
        op.runtleap(work_dir=directory + "/", mol_file='mol.pdb', tleaptemp=settings.tleap_template, tleapin="leap.in")

        op.runAmberMPI(work_dir=directory + "/", np=settings.mpi_procs, amberin=settings.amber_params)
        
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
        
        finalEnergy = op.parseAmberEnergy(directory + "/amber.out")
        
        individuals[i].setFitness(fitness_index, finalEnergy)
        
        if (finalEnergy == 999):
            continue
        
        obConversion = openbabel.OBConversion() 
        obConversion.SetInFormat("pdb")

        molec = openbabel.OBMol()
        obConversion.ReadFile(molec, directory + '/min_struct.pdb')
        individuals[i].mol = molec
    

        print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))
        
def amber_energy(settings, individuals, fitness_index, pop_start=0):

    # os.environ["AMBERHOME"]="/"+settings.amber_path

    # if (settings.verbose):
        # print("IND NUM, IND FITNESS")
        
    for i in xrange(pop_start, len(individuals)):
    
        # Updating the OBmol of each individual with current dihedrals
        if (settings.composition_optimization and not individuals[i].init):
            individuals[i].applyComposition(settings)

        if (settings.dihedral_optimization and not individuals[i].init):
            individuals[i].applyPhiPsiDihedrals()

        if (settings.basilisk_and_sidechains and not individuals[i].init):
            individuals[i].apply_chi_dihedrals()

            
        # Writing the molecule for further Amber input in a dedicated folder
        # directory = settings.output_path + "/amber_run/ind_"+str(i)
        directory = settings.output_path + "/amber_run"
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Seperate if non-standard residues are involved
        if (not settings.use_gaussian):
        
            if (os.path.exists(directory + "/in_temp_mol.pdb")):
                os.remove(directory + "/in_temp_mol.pdb")

            obconv = openbabel.OBConversion()
            obconv.SetOutFormat("pdb")
            obconv.WriteFile(individuals[i].mol, directory + "/in_temp_mol.pdb")

            '''
            # Preprocessing pdb file if needed
            # Creating a more suitable pdb file for tleap to avoid further issues
            # pdb4amber analyses PDB files and cleans them for further usage, especially with the LeaP programs of Amber
            if (os.path.exists(directory+"/pdb4amber_temp_mol.pdb")):
                os.remove(directory+"/pdb4amber_temp_mol.pdb")

            foutanderr = open(directory+'/pdb4amber_temp_mol.log','w')
            try:
                proc = subprocess.Popen(["pdb4amber", "-i", directory+"/in_temp_mol.pdb", "-o", directory+"/pdb4amber_temp_mol.pdb", "--dry", "--nohyd"], stdout=foutanderr, stderr=subprocess.STDOUT) #"--nohyd", 
                proc.wait()
                foutanderr.close()
            except IOError as e:
                sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
            except subprocess.CalledProcessError as e:
                sys.exit("pdb4amber failed, returned code %d (check '"+directory+"/pdb4amber_temp_mol.log')" % (e.returncode))
            except OSError as e:
                sys.exit("failed to execute pdb4amber: %s" % (str(e)))
            '''
            '''    
            #(Re-)adding the hydrogen atoms if "--nohyd" flag used to remove hydrogens with pdb4amber previously, actually reduce can also be used in pdb4amber with --reduce, this is used as tleap is very strict on hydrogens' names
            if (os.path.exists(directory+"/reduce_temp_mol.pdb")):
                os.remove(directory+"/reduce_temp_mol.pdb")

            fout = open(directory+'/reduce_temp_mol.pdb','w')
            ferr = open(directory+'/reduce_temp_mol.log', 'w')
            try:
                proc = subprocess.Popen(["reduce", "-nuclear", "-NOCon", directory+"/in_temp_mol.pdb"], stdout=fout, stderr=ferr)# "-build", "-nuclear"
                proc.wait()
                fout.close()
                ferr.close()
            except IOError as e:
                sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
            except subprocess.CalledProcessError as e:
                sys.exit("reduce failed, returned code %d (check '"+directory+"/reduce_temp_mol.log')" % (e.returncode))
            except OSError as e:
                sys.exit("failed to execute reduce: %s" % (str(e)))
            '''

            # Preprocessing the molecule with antechamber
            if (os.path.exists(directory + "/antechamber_temp_mol.mol2")):
                os.remove(directory + "/antechamber_temp_mol.mol2")

            foutanderr = open(directory + '/antechamber_temp_mol.log', 'w')
            try:
                proc = subprocess.Popen(["antechamber", "-i", directory + "/in_temp_mol.pdb", "-fi", "pdb", "-o", directory + "/antechamber_temp_mol.mol2", "-fo", "mol2", "-pf", "y"], stdout=foutanderr, stderr=subprocess.STDOUT)  # skips optimization before AM1-BCC charge computation / "-ek", "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0"/ "-rn", "MOL"/"-nc", "0 / -ek,  ndiis_attempts=700"
                proc.wait()
                foutanderr.close()  # "-c", "bcc", "-ek", "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0"
            except IOError as e:
                sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
            except subprocess.CalledProcessError as e:
                sys.exit("antechamber failed, returned code %d (check '" + directory + "/pdb4amber_temp_mol.log')" % (e.returncode))
            except OSError as e:
                sys.exit("failed to execute antechamber: %s" % (str(e)))
                
            tleap_input_mol = "in_temp_mol.pdb"
            
        
        elif (settings.use_gaussian):
        
            if (os.path.exists(directory + "/in_temp_mol.com")):
                os.remove(directory + "/in_temp_mol.com")

            obconv = openbabel.OBConversion()
            obconv.SetOutFormat("com")
            obconv.WriteFile(individuals[i].mol, directory + "/in_temp_mol.com")
            '''
            # Converting the pdb molecule to Gaussian input file, often buggy
            foutanderr = open(directory+'/gaussian_temp_mol.log','w')
            try:
            proc = subprocess.Popen(["newzmat", "-ipdb", "-ocart", "-tcart", "-sMMC2", directory+"/in_temp_mol.pdb", directory+"/in_temp_mol.com", settings.output_path + "/gaussian_template.com"], stdout=foutanderr, stderr=subprocess.STDOUT)
               
            '''

            # For non-standard residues or use of Gaussian for charge computation
            
            # Optimizing the structure first if needed
            # Preparing Gaussian input file
            '''
            proccp = subprocess.Popen(['cp', settings.output_path+"/opt_gaussian_template.gin", directory+"/opt_gaussian_temp_mol.gin"], stdout=subprocess.PIPE)
            out, err = proccp.communicate()

            f = open(directory+"/in_temp_mol.com")
            lines = f.readlines()
            with open(directory+"/opt_gaussian_temp_mol.gin", "a") as myfile:
                myfile.writelines(lines[5:])
            f.close()
            myfile.close()

            # Running Gaussian to optimize geometry
            proc = subprocess.Popen(["g16",  directory+"/opt_gaussian_temp_mol.gin"], stdout=subprocess.PIPE)
            out, err = proc.communicate()
            '''

            # Running single-point Gaussian to compute charges for non-standard residues from optimized geometry 
            proccp = subprocess.Popen(['cp', settings.output_path + "/sp_gaussian_template.gin", directory + "/sp_gaussian_temp_mol.gin"], stdout=subprocess.PIPE)
            out, err = proccp.communicate()
            
            f = open(directory + "/in_temp_mol.com")
            lines = f.readlines()
            with open(directory + "/sp_gaussian_temp_mol.gin", "a") as myfile:
                myfile.writelines(lines[5:])
            f.close()
            myfile.close()
            
            proc = subprocess.Popen(["g16", directory + "/sp_gaussian_temp_mol.gin"], stdout=subprocess.PIPE)
            out, err = proc.communicate()


            # Preprocessing the molecule with antechamber, especially for non-standard residues and partial charges
            if (os.path.exists(directory + "/antechamber_temp_mol.mol2")):
                os.remove(directory + "/antechamber_temp_mol.mol2")

            foutanderr = open(directory + '/antechamber_temp_mol.log', 'w')
            try:
                proc = subprocess.Popen(["antechamber", "-i", directory + "/sp_gaussian_temp_mol.log", "-fi", "gout", "-o", directory + "/antechamber_temp_mol.mol2", "-fo", "mol2", "-c", "resp"], stdout=foutanderr, stderr=subprocess.STDOUT)
                proc.wait()
                foutanderr.close()
            except IOError as e:
                sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
            except subprocess.CalledProcessError as e:
                sys.exit("antechamber failed, returned code %d (check '" + directory + "/antechamber_temp_mol.log')" % (e.returncode))
            except OSError as e:
                sys.exit("failed to execute antechamber: %s" % (str(e)))

                
            tleap_input_mol = "antechamber_temp_mol.mol2"


        # Checking missing parameters with parmchk and generate the frcmod with missing params (only once since the frcmod file should remain the same for all energy computations)
        if (os.path.exists(directory + "/parmchk_temp_mol.frcmod")):
            # os.remove(directory+"/parmchk_temp_mol.frcmod")
            pass
        
        else:
            
            foutanderr = open(directory + '/parmchk_temp_mol.log', 'w')
            try:
                proc = subprocess.Popen(["parmchk2", "-i", directory + "/antechamber_temp_mol.mol2", "-f", "mol2", "-o", directory + "/parmchk_temp_mol.frcmod"], stdout=foutanderr, stderr=subprocess.STDOUT)
                proc.wait()
                foutanderr.close()
            except IOError as e:
                sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
            except subprocess.CalledProcessError as e:
                sys.exit("parmchk2 failed, returned code %d (check '" + directory + "/parmchk_temp_mol.log')" % (e.returncode))
            except OSError as e:
                sys.exit("failed to execute parmchk2: %s" % (str(e)))
        '''    
        procmv = subprocess.Popen(['cp', directory+"/parmchk_temp_mol.frcmod", directory+"/parmchk_{}.frcmod".format(i)], stdout=subprocess.PIPE) # see also 'logfile' if pmemd.MPI is used 
        out, err = procmv.communicate()
        '''

          
        # Preparing input files with tleap, use pdb4amber/reduce/in_temp_mol.pdb depending on the pdb file relevance
        op.runtleap(work_dir=directory + "/", mol_file=tleap_input_mol, tleaptemp=settings.tleap_template, tleapin="leap.in")

        
        # Computing the energy of the molecule with Amber/Pmemd
        op.runAmberMPI(work_dir=directory + "/", np=settings.mpi_procs, amberin=settings.amber_params)
        

        # Generate output minimized molecule with amber in a pdb file
        fout = open(directory + '/amber_out_min_struct.pdb', 'w')
        ferr = open(directory + '/amber_out_min_struct.log', 'w')
        try:
            proc = subprocess.Popen(["ambpdb", "-p", directory + "/leap_temp_mol.prmtop", "-c", directory + "/amber.rst"], stdout=fout, stderr=ferr)
            # out, err = proc.communicate(input=directory+"/amber.rst")
            # fout.write(out)
            # ferr.write(err)
            proc.wait()
            fout.close()
            ferr.close()
        except IOError as e:
            sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
        except subprocess.CalledProcessError as e:
            sys.exit("convert rst to pdb failed, returned code %d (check '" + directory + "/amber_out_min_struct.log')" % (e.returncode))
        except OSError as e:
            sys.exit("failed to convert rst to pdb: %s" % (str(e)))

        ''' If simulated in a solvent, to remove the solvent    
        # use pdb4amber and reduce again to remove water from the output molecule but keep the hydrogens if needed
        if (os.path.exists(directory+"/out_min_mol.pdb")):
            os.remove(directory+"/out_min_mol.pdb")
            
        foutanderr = open(directory+'/out_min_mol.log','w')
        try:
            proc = subprocess.Popen(["pdb4amber", "-i", directory+"/amber_out_min_struct.pdb", "-o", directory+"/out_min_mol.pdb", "--dry"], stdout=foutanderr, stderr=subprocess.STDOUT)
            proc.wait()
            foutanderr.close()
        except IOError as e:
            sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
        except subprocess.CalledProcessError as e:
            sys.exit("pdb4amber failed, returned code %d (check '"+directory+"/out_min_mol.log')" % (e.returncode))
        except OSError as e:
            sys.exit("failed to execute pdb4amber: %s" % (str(e)))
        '''    

        # Moving some output files generated by amber until now
        procmv = subprocess.Popen(['mv', 'leap.log', 'mdinfo', directory], stdout=subprocess.PIPE)  # see also 'logfile' if pmemd.MPI is used 
        out, err = procmv.communicate()



        # Random reseeding within population if bad configurations occur in amber (explosion of system during minimization), otherwise normal update of molecule and energy
        error_detected = False
        f = open(directory + "/amber.out", "r")
        filelines = f.readlines()
        
        for line in filelines:
            if ('SANDER BOMB') in line:
                error_detected = True
                
            elif ('NaN') in line:
                error_detected = True
                      
            # elif ('VDWAALS = ******') in line:
                # error_detected = True
                
        f.close()

        if (error_detected):
            '''
            from random import randint
            new_index = randint(0, settings.population_size-1)
            if new_index == i:
                new_index = new_index - 1

            print('ILL INDIVIDUAL, re-seeding indiv {} randomly replaced by current indiv {}'.format(i, new_index))

            individuals[i].mol = individuals[new_index].mol
            individuals[i].setFitness(fitness_index, individuals[new_index].fitnesses)
            '''
            print('ILL INDIVIDUAL, re-seeding indiv {} according to initilization settings'.format(i))

            if (settings.mc_generate_dihedrals):
                # "Generating Initial Population dihedrals via MC"
                generators.MonteCarloDihedralGenerator(settings, [individuals[i]])

            else:
                # "Generating Initial Population dihedrals from uniform distribution"
                generators.UniformDihedralGenerator(settings, [individuals[i]])
                
            if (settings.basilisk_and_sidechains):
                # "Backbone-dependent side chains generated by Basilisk"
                generators.BasiliskSideChainDihedralsGenerator(settings, [individuals[i]])
            
                
            amber_energy(settings, individuals, fitness_index, pop_start=i)
            break

        
        elif (not error_detected):

            # Parse the energy in output file and assign it to the individual's fitness 
            finalEnergy = op.parseAmberEnergy(directory + "/amber.out")

            # Detection of energy explosion, defined threshold, to allow OpenBabel to load the output file with connectivity based on distances between atoms 
            if finalEnergy > 500:
                
                print('ENERGY EXPLOSION, re-seeding indiv {} according to initilization settings'.format(i))

                if (settings.mc_generate_dihedrals):
                    # "Generating Initial Population dihedrals via MC"
                    generators.MonteCarloDihedralGenerator(settings, [individuals[i]])

                else:
                    # "Generating Initial Population dihedrals from uniform distribution"
                    generators.UniformDihedralGenerator(settings, [individuals[i]])

                if (settings.basilisk_and_sidechains):
                    # "Backbone-dependent side chains generated by Basilisk"
                    generators.BasiliskSideChainDihedralsGenerator(settings, [individuals[i]])


                amber_energy(settings, individuals, fitness_index, pop_start=i)
                break
        
            else:
        
                # Assign the output molecule to the individual's molecule
                # Read pdb file
                path = str((directory + '/amber_out_min_struct.pdb').strip('\''))
                obConversion = openbabel.OBConversion()  # openbabel parse
                obConversion.SetInFormat(path.split(".")[1])

                molec = openbabel.OBMol()
                if (not obConversion.ReadFile(molec, path)):
                    print("Problem reading " + path)
                    sys.exit(0)

                individuals[i].mol = molec
                individuals[i].setFitness(fitness_index, finalEnergy)

                print('ind {}, fitness {}'.format(i, individuals[i].fitnesses))

                # Update torsion angles values of the individual
                phi_psi_dihedrals = individuals[i].getPhiPsiDihedrals()

                for j in xrange (0, len(phi_psi_dihedrals)):
                    individuals[i].phi_dihedrals[j] = phi_psi_dihedrals[j][0]
                    individuals[i].psi_dihedrals[j] = phi_psi_dihedrals[j][1]

                if (settings.basilisk_and_sidechains):
                    chi_dihedrals = individuals[i].get_chi_dihedrals_per_res()
                    individuals[i].chi_angles = chi_dihedrals
