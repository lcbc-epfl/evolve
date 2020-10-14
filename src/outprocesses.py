'''
All processes that operate outside of the EVOLVE main code like simulation codes should be included in this.

.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu
.. codeauthor:: Justin Villard
'''

from __future__ import division
from __future__ import print_function

from builtins import str
from builtins import range
from past.utils import old_div
import sys

import os
import subprocess
import numpy as np
import openbabel
import time


def runtleap(work_dir="", mol_file="mol.pdb", tleaptemp="tleap_template.in", tleapin="leap.in", inpcrd_out="mol.inpcrd", prmtop_out="mol.prmtop"):
    '''

    runs tleap to process a pdb file into a valid amber prmtop file.
    This procedure needs to be verified before runinng evolve such that the topology building process does not fail.
    If the setting for the resname is not set the mol files need to be reloaded.

    Parameters
    ----------

    work_dir : string
        working directory
    mol_file : string
        path to PDB file
    tleaptemp
    tleapin
    inpcrd_out
    prmtop_out

    Returns
    -------

    '''
    
    if (os.path.exists(work_dir + inpcrd_out)):
        os.remove(work_dir + inpcrd_out)
    if (os.path.exists(work_dir + prmtop_out)):
        os.remove(work_dir + prmtop_out)
        
    f = open(tleaptemp, "r")
    g = f.readlines()
    filename = (work_dir + mol_file).split('.')[0]
    
    for i, line in enumerate(g):
        if ("mol" and "=") in line and ("mol2" in mol_file):
            g[i] = "mol = loadmol2 " + work_dir + mol_file + "\n"
        elif ("mol" and "=") in line and ("pdb" in mol_file):
            g[i] = "mol = loadpdb " + work_dir + mol_file + "\n"
        elif ("saveoff") in line:
            g[i] = "saveoff mol " + filename + ".lib\n"
        elif "saveamberparm" in line:
            g[i] = "saveamberparm mol " + work_dir + prmtop_out + " " + work_dir + inpcrd_out + "\n"
    
    ff = open(work_dir + tleapin, 'w')
    for line in g:
        ff.write(line)
    ff.close()
    f.close()

    # proc = subprocess.Popen([tleap_path+"/bin/tleap", "-f", work_dir+tleapin], stdout=subprocess.PIPE)
    # out, err = proc.communicate()

    if (work_dir == ""):
        foutanderr = open('leap.log', 'w')
    else:
        foutanderr = open(work_dir + '/leap.log', 'w')
    try:
        proc = subprocess.Popen(["tleap", "-f", work_dir + tleapin], stdout=foutanderr, stderr=subprocess.STDOUT)
        proc.wait()
        foutanderr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("tleap failed, returned code %d (check '" + work_dir + "/leap.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute tleap: %s" % (str(e)))


def runtleapMultiInd(work_dir="", leap_inputfile="tleap_template.in", n_frames=10):
    '''

    Takes the pdb file written to disk after openbabel mutation,
    then uses pdb4amber to process it and add TER statements for tleap.
    Then tleap is run on all frames of the individual.

    If the option `use_res_type=True` is not used, then you need to include all necessary *.off files for the a99SB forcefield in the share directory
    If you want to track which rotamer was placed this might be an option. If you use a forcefield different than a99SB you have to use `use_res_type=True`

    Parameters
    ----------
    work_dir: str
        working directory
    leap_inputfile: str
        path to inputfile in config file
    n_frames: int
        number of frames

    '''
    logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
    logfile.write("1. Running parameter generation for " + str(n_frames) + " frames\n")
    logfile.write("1.1 Converting files from openbabel to amber format\n")
    processes = []
    try:
        for frame in range(0, n_frames):
            os.chdir(work_dir + "/frame" + str(frame))
            f = open('mutated_ter.pdb', "w")
            FNULL = open(os.devnull, 'w')
            pPDB4Amber = subprocess.Popen(["pdb4amber", "mutated.pdb", "--nohyd"], universal_newlines=True, stdout=f,
                                          stderr=FNULL)
            processes.append((pPDB4Amber, f, FNULL))
            # f.close()  # In order to avoid I/O errrors  with low ulimit it might be necessary to use subprocess call. Then files should be immediateluy closed
            # FNULL.close()
            os.chdir("../../")
            pass
        for p, f, FNULL in processes:
            p.wait()
            FNULL.close()
            f.close()
            pass
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("pdb4amber failed processing mutated.pdb, returned code %d " % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute pdb4amber: %s" % (str(e)))
    pass
    # then run leap in every directory, this will generate all the necessary parameter files
    logfile.write("1.2 Running tleap\n")
    processes = []
    try:
        for frame in range(0, n_frames):
            os.chdir(work_dir + "/frame" + str(frame))
            foutanderr = open('leap.log', 'w')
            pLeap = subprocess.Popen(["tleap", "-f", leap_inputfile], stdout=foutanderr, stderr=subprocess.STDOUT)
            processes.append((pLeap, foutanderr))
            os.chdir("../../")
            pass
        for p, f in processes:
            p.wait()
            f.close()
            pass
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("tleap failed, returned code %d (check '" + work_dir + "/leap.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute tleap: %s" % (str(e)))
    pass
    logfile.close()


def runPMEMDMultiInd(work_dir="", np=20, amberin='', n_frames=10, compute_cluster=False, ntasks=60, nodes=6, queuename="debug"):
    '''

    minimization of multiple frames in parallel for mmpbsa_multi evuluator
    All files are put in work_dir/frame* folders
    We're using the multipmemd version of pmemd.MPI. For this we generate a groupfile in the work dir and then call pmemd
    Cannot use pmemd.cuda because of GPU precision model. With unreasonable geometry (atoms very close together which can be the case if swapping a large sidechain for a small sidechain with the GA) forces can be truncated or, in the worst case, overflow the fixed precision presentation.

    In order to maximize the performance the multipmemd module is used.

    #ToDO
    Gromacs should be possible here as well using parmed

    Parameters
    ----------
    work_dir: str
        working directory
    np: int
        number of processor cores
    amberin: : str
        /path/to/amberinput file
    n_frames: int
        number of frames
    compute_cluster: bool
        Determines whether a submit script shall be created and submitted using sbatch
    ntasks: int
        number of tasks for slurm
    nodes: int
        number of nodes for slumr
    queuename: str
        name of the queue

    '''

    logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
    logfile.write("Ran minimization using following command: ")
    f = open(work_dir + "/groupfile", "w")
    f.write('')
    f.close()
    f = open(work_dir + "/groupfile", "a")
    for frame in range(0, n_frames):
        path = work_dir + "/frame" + str(frame) + "/"
        f.write(
            '-O -i ' + amberin + ' -o ' + path + 'min.out -p ' + path + 'GB1AZI_solv.prmtop -c ' + path + 'GB1AZI_solv.inpcrd -r ' + path + 'min.rst -suffix frame' + str(
                frame) + ' -l ' + path + 'logfile -inf ' + path + 'mdinfo \n')
        # f.write('-O -i '+amberin+' -o '+path+'min.out -p '+path+'GB1AZI_solv.prmtop -c '+path+'GB1AZI_solv.inpcrd -r '+path+'min.rst -suffix frame'+str(frame)+"\n") #for use in sander
        pass
    f.close()
    individualStart = time.time()
    # If we calculate on a compute cluster make runscript and run command
    # potentially it should be changed from sbatch to srun to have the complete GA in a slurm script to not block one core on the login node
    if (compute_cluster):
        if (os.path.isfile('./run_pmemd.sh') == False):
            pmemdscript = "#!/bin/bash\n#SBATCH -J GA-min\n#SBATCH --partition={}\n#SBATCH --time=00:05:00\n#SBATCH --nodes={}\n#SBATCH --ntasks-per-node={}\n#SBATCH --mem=10000\n#SBATCH --output=slurm/job.out\n#SBATCH --error=slurm/error.err\n### Load Modules ###\nmodule load amber/16/intel-17.0.4\nsrun pmemd.MPI -ng {} -groupfile {} \n".format(
                str(queuename), str(nodes), str(ntasks), str(n_frames), work_dir + "/groupfile")
            pmemdscript_f = open("run_pmemd.sh", "w")
            pmemdscript_f.write(pmemdscript)
            pmemdscript_f.close()
            pass
        # execute with sbatch and wait for command to finish
        pMin = subprocess.Popen('sbatch -W run_pmemd.sh', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = pMin.communicate()
        pMin.wait()
        if (pMin.returncode != 0 and pMin.returncode != 174):
            sys.exit("failed to execute sbatch and pmemd.MPI. Only slurm is supported\n: %d, %s, %s" % (
            pMin.returncode, output, error))
        if pMin.returncode == 174:
            print("Segmentation fault occured, skipping this individual")
            return 0
            pass
        pass
    else:
        logfile.write('mpirun -np ' + str(np) + ' pmemd.MPI -ng ' + str(n_frames) + ' -groupfile ' + work_dir + '/groupfile\n')
        pMin = subprocess.Popen('mpirun -np ' + str(np) + ' pmemd.MPI -ng ' + str(n_frames) + ' -groupfile ' + work_dir + '/groupfile',stdout=logfile, shell=True)
        output, error = pMin.communicate()
        pMin.wait()
        if (pMin.returncode != 0 and pMin.returncode != 174):
            sys.exit("failed to execute sbatch and pmemd.MPI. Only slurm is supported\n: %d, %s, %s" % (
            pMin.returncode, output, error))
        if pMin.returncode == 174:
            print("Segmentation fault occured, skipping this individual")
            return 0
            pass
    individualEnd = time.time()
    logfile.write("Runtime of pmemd " + str(individualEnd - individualStart) + "sec." + "\n")
    # after pmemd has finished we convert the rst file to a pdb file
    logfile.write("Converting rst files to pdb format\n")
    processes = []
    try:
        for frame in range(0, n_frames):
            path = work_dir + '/frame' + str(frame) + '/'
            f = open(path + 'min.pdb', "w")
            pPDBConv = subprocess.Popen(['ambpdb', '-p', path + 'GB1AZI_solv.prmtop', '-c', path + 'min.rst.frame' + str(frame)],universal_newlines=True, stdout=f)
            processes.append((pPDBConv, f))
            pass
        for p, f in processes:
            p.wait()
            f.seek(0)
            f.close()
            pass
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("ambpdb failed, returned code %d" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute ambpdb: %s" % (str(e)))
    pass
    # then we remove the solvent waters from the minimization, This is necessary as the different frames can have different amounts of solvent and therefore the prmtop files are different, which cannot be the case for the MMPBSA calculation.
    logfile.write("Desolvating complex using pdb4amber -d and removing counter ions\n")
    processes = []
    individualStart = time.time()
    try:
        for frame in range(0, n_frames):
            os.chdir(work_dir + "/frame" + str(frame))
            f = open("min_nosolv.pdb", "w")
            FNULL = open(os.devnull, 'w')
            pSolvStrip = subprocess.Popen(["pdb4amber", "-s", ":Cl-,Na+", "-d", "min.pdb"], universal_newlines=True,stdout=f, stderr=FNULL)
            # f.close()  # In order to avoid I/O errrors use subprocess call.
            # FNULL.close()
            processes.append((pSolvStrip, f, FNULL))
            os.chdir("../../")
            pass

        for p, f, fnull in processes:
            p.wait()
            f.close()
            fnull.close()
            pass
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("pdb4amber -d failed, returned code %d " % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute pdb4amber -d: %s" % (str(e)))

    individualEnd = time.time()
    logfile.write("Runtime of Desolvation: " + str(individualEnd - individualStart) + " sec." + "\n")
    logfile.close()
    return 1
    pass


def runMMPBSAMultiInd(work_dir="", mmpbsa_inputfile="mmpbsa.in", n_frames=10, compute_cluster=False, ntasks=60, nodes=6,queuename="debug", energy_evaluator="gb"):
    '''

    MMPBSA calculation with mmpbsa_multi evaluator. Takes all minimized files for the frames of the Individual and runs an MMPBSA calclation using MMPBSA.py.MPI
    Extracts the MMPBSA energy and returns it.
    Stores the energy from the implicit solvent calculation in the `work_dir/evaluators_stability_multi.log`

    Parameters
    ----------
    work_dir: str
        working directory
    mmpbsa_inputfile: str
        /path/to/mmpbsa input file
    n_frames: int
        number of frames
    compute_cluster: bool
        Determines whether a submit script shall be created and submitted using sbatch
    ntasks: int
        number of tasks for slurm
    nodes: int
        number of nodes for slumr
    queuename: str
        name of the queue

    Returns
    -------
    binding_affinity: float
        value extracted from MMPBSA.log file
    '''
    framesstring = ""
    logfile = open(work_dir + "/evaluator_mmpbsa_multi.log", "a")
    # construct framestring with paths to all frames
    for frame in range(0, n_frames):
        framesstring += "./frame" + str(frame) + "/min_nosolv.pdb "
    pass
    # If we calculate on a compute cluster make runscript and run command
    # potentially it should be changed from sbatch to srun to have the complete GA in a slurm script to not block one core on the login node
    os.chdir(work_dir)
    if (compute_cluster):
        if (os.path.isfile('./run_mmpbsa.sh') == False):
            mmpbsascript = """#!/bin/bash
#SBATCH -J GA-MMPBSA
#SBATCH --partition=%s
#SBATCH --time=00:05:00
#SBATCH --nodes=%s
#SBATCH --ntasks-per-node=%s
#SBATCH --mem=10000
#SBATCH --output=slurm/job.out
#SBATCH --error=slurm/error.err
module load amber/16/intel-17.0.4
srun -n %s MMPBSA.py.MPI -O -i inp/mmpbsa.in -o work_dir/mmpbsa.log -cp work_dir/frame0/dry_complex.prmtop -rp work_dir/frame0/receptor.prmtop -lp work_dir/frame0/ligand.prmtop -y %s""" % (
            str(queuename), str(nodes), str(ntasks), str(n_frames), framesstring)
            mmpbsascript_f = open("run_mmpbsa.sh", "w")
            mmpbsascript_f.write(mmpbsascript)
            mmpbsascript_f.close()
            pass
        logfile.write("Submitting batch job and wait for it to finish\n")
        pMMPBSA = subprocess.Popen("sbatch -W run_mmpbsa.sh", stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   shell=True)
        output, error = pMMPBSA.communicate()
        pMMPBSA.wait()
        if pMMPBSA.returncode != 0:
            sys.exit("failed to execute sbatch and run MMPBSA calculation.\n: %d, %s, %s" % (
            pMMPBSA.returncode, output, error))
        pass
        pass
    else:
        logfile.write("mpirun -np " + str(
            n_frames) + " MMPBSA.py.MPI -O -i " + mmpbsa_inputfile + " -o /mmpbsa.log -cp ./frame0/dry_complex.prmtop -rp ./frame0/receptor.prmtop -lp ./frame0/ligand.prmtop -y " + framesstring + "\n\n")
        logfile.close()  # need to open and reclose to have correct order of writing
        logfile = open("./evaluator_mmpbsa_multi.log", "a")
        pMMPBSA = subprocess.Popen("mpirun -np " + str(
            n_frames) + " MMPBSA.py.MPI -O -i " + mmpbsa_inputfile + " -o ./mmpbsa.log -cp ./frame0/dry_complex.prmtop -rp ./frame0/receptor.prmtop -lp ./frame0/ligand.prmtop -y " + framesstring,
                                   stdin=None, stdout=logfile, stderr=None, shell=True)
        output, error = pMMPBSA.communicate()
        pMMPBSA.wait()
        if pMMPBSA.returncode != 0:
            sys.exit("failed to execute MMPBSA.MPI : %d, %s, %s" % (pMMPBSA.returncode, output, error))
        pass
    os.chdir("../")
    # Grab output from mmpbsa calculation
    fp = open(work_dir + "/mmpbsa.log")
    if energy_evaluator == "gb":
        linenumber = 87 + n_frames  # GB
        pass
    else:
        linenumber = 150 + n_frames  # GB+PB calculation line number
        pass
    for i, line in enumerate(fp):
        if i == int(linenumber):  # set to 96 for gb or 159 to pb calculation
            binding_affinity_GB = line[25:36]
            break
    fp.close()
    logfile.close()
    return float(binding_affinity_GB)
    pass


def writeFrontLog(path="", front='', iteration=0):
    '''

    write front log as semicolon seperated values with fitnesses seperated by comma

    Parameters
    ----------
    path: str
        path to the current working directory, not work_dir/
    front: int
        id of the front
    iteration: int
        iteration counter of the GA
    '''
    filename = path + '/first_fronts.log'

    print("Writing first front to", filename)

    if iteration > 0:
        append_write = 'a'  # append if already exists
    else:
        append_write = 'w'  # make a new file if not
    string = ''
    for item in front:
        string += str(item[0]) + ',' + str(item[1]) + ';'
    frontsfile = open(filename, append_write)
    frontsfile.write(string + '\n')
    frontsfile.close()
    pass


def writeEnergyLog(path="", energies=None, iteration=0):
    '''

    write log of energies for all min_mols as comma seperated values 

    The format is energy wt, correction wt, energy_mut, correction_mut

    Parameters
    ----------
    path: str
        path to the current working directory, not work_dir/
    energies: list
        (energy wt, correction wt, energy_mut, correction_mut)
    iteration: int
        iteration counter of the GA
    '''
    filename = path + '/energy.log'

    print("Writing energies", filename)

    if iteration > 0:
        append_write = 'a'  # append if already exists
    else:
        append_write = 'w'  # make a new file if not
    string = ''
    for item in energies.keys():
        string += str(energies[item]) + ','
    frontsfile = open(filename, append_write)
    frontsfile.write(string + '\n')
    frontsfile.close()
    pass


def writeFrontFile(path="", file='', fitness='', iteration=0, frontid=0):
    '''

    write pdb files for each molecule in the first pareto front

    Parameters
    ----------
    path: str
        path to the current working directory, not work_dir/
    file: Object
        OBMol to write to disk
    fitness: list
        list of fitness values
    iteration: int
        iteration counter of the GA
    frontid: int
        id of the individual in the front to write to disk
    '''
    if not os.path.exists(path + "/intermediate_molecules"):
        os.makedirs(path + "/intermediate_molecules")
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdb")
    file.SetTitle("fitness:" + fitness)
    obConversion.WriteFile(file, path + '/intermediate_molecules/' + "min_iter_" + str(iteration + 1) + "_front_" + str(
        frontid) + ".pdb")

    pass


def runAmberMPI(work_dir="", np=16, amberin="amber_minimization.in", prmtop="mol.prmtop", inpcrd="mol.inpcrd", amberout="amber.out", restart="amber.rst"):
    '''

    Amber MPI documentation

    %TODO
    Why is this there? Should go into one block with pmemd

    Parameters
    ----------
    work_dir: str
        path to calculations where amber was run
    np: int
        number of cores, set in the input file
    amberin: str
        file to amber input
    prmtop: str
        name of the amber parameter file
    inpcrd: str
        name of the input coordinates
    amberout: str
        name of amber out put file
    restart: str
        name of restart file

    Returns
    -------

    '''

    # Actually already made by amber with flag -O
    #if (os.path.exists(work_dir + amberout)):
    #    os.remove(work_dir + amberout)
    #if (os.path.exists(work_dir + restart)):
    #    os.remove(work_dir + restart)
    #
    # proc = subprocess.Popen(["mpirun", "-np", str(np), "sander.MPI", "-O", "-i", amberin, "-o", work_dir+amberout, "-c", work_dir+inpcrd, "-p", work_dir+prmtop, "-r", work_dir+restart], stdout=subprocess.PIPE)
    # can use pmemd if not in vaccum, more efficient parallelization but apparently doesn't work in vaccum
    # out, err = proc.communicate()
    
    if (work_dir == ""):
        foutanderr = open('amber.log', 'w')
    else:
        foutanderr = open(work_dir + '/amber.log', 'w')

    try:
        if int(np) > 1:
            proc = subprocess.Popen(["mpirun", "-np", str(np), "sander.MPI", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        else:
            proc = subprocess.Popen(["sander", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        proc.wait()
        foutanderr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("sander failed, returned code %d (check '" + work_dir + "/amber.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute sander.MPI: %s" % (str(e)))

    
def runPMEMD(work_dir="", np=20, amberin="amber_minimization.in", prmtop="mol.prmtop", inpcrd="mol.inpcrd", amberout="amber.out", restart="amber.rst", use_cuda=False):
    '''

    run PMED using the specified number of cores .
    use cuda should be set to false for all minimizations/ early heating!
    Otherwise forces might be truncated


    Parameters
    ----------
    work_dir: str
        path to calculations where amber was run
    np: int
        number of cores, set in the input file
    amberin: str
        file to amber input
    prmtop: str
        name of the amber parameter file
    inpcrd: str
        name of the input coordinates
    amberout: str
        name of amber out put file
    restart: str
        name of restart file
    use_cuda: boolean
        determine if cuda should be used

    Returns
    -------
    nothing, exits if amber has problems
    '''
    
    if (work_dir == ""):
        foutanderr = open('amber.log', 'w')
    else:
        foutanderr = open(work_dir + '/amber.log', 'w')

    try:
        if use_cuda:
            proc = subprocess.Popen(["pmemd.cuda", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        elif int(np) > 1:
            proc = subprocess.Popen(["mpirun", "-np", str(np), "pmemd.MPI", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        else:
            proc = subprocess.Popen(["pmemd", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        proc.wait()
        foutanderr.close()
    
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("sander failed, returned code %d (check '" + work_dir + "/amber.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute PMEMD.MPI: %s" % (str(e)))
    

def runMMPBSA(work_dir="", mmpbsa_in="mmpbsa.in", prmtop="mol.prmtop", inpcrd="mol.inpcrd", final_results="mmpbsa.dat", decomp_results="mmpbsa_decomp.dat"):
    '''

    runs a mmpbsa calculation using the files in the workdir.
    File names are defined in the input file
    Runs a decomposition analysis too.

    Used with evaluator decomp

    Parameters
    ----------
    work_dir: str
        path to calculations where amber was run
    mmpbsa_in: str
        name of the mmpbsa file
    prmtop: str
        name of the amber parameter file
    inpcrd: str
        name of the input coordinates
    final_results:str
        name of the mmpbsa dat file
    decomp_results: str
        name of the file for the decompositon output

    Returns
    -------
    nothing
    '''
    
    if (work_dir == ""):
        foutanderr = open('mmpbsa.log', 'w')
    else:
        foutanderr = open(work_dir + '/mmpbsa.log', 'w')
    
    try:
        proc = subprocess.Popen(["MMPBSA.py", "-O", "-i", mmpbsa_in, "-o", work_dir + final_results, "-do", work_dir + decomp_results, "-cp", work_dir + prmtop, "-y", work_dir + inpcrd], stdout=foutanderr, stderr=subprocess.STDOUT)
        proc.wait()
        foutanderr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("MMPBSA failed, returned code %d (check '" + work_dir + "/mmpbsa.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute MMPBSA.py: %s" % (str(e)))


def parseMMPBSA(mmpbsa_file_path):
    '''

    read MMPBSA energy


    Parameters
    ----------
    mmpbsa_file_path

    Returns
    -------

    '''
    import os
    
    if not os.path.isfile(mmpbsa_file_path):
        return 999.0
    
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Energy Component" in lines[i]):
            start = i
            break
        
    if (start == -1):
        return
    
    return float(lines[start + 15].split()[1])

    
def parseMMBPSA_total_decomp(mmpbsa_file_path, num_res):
    '''

    parse MMPBSA total decompostion

    #TODO
    Nick

    Parameters
    ----------
    mmpbsa_file_path
    num_res

    Returns
    -------

    '''
    import os
    
    if not os.path.isfile(mmpbsa_file_path):

        return None
    
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Total Energy Decomposition" in lines[i]):
            start = i
            break
        
    if (start == -1):
        return
    
    curr = i + 3
    end = i + 3 + (num_res)
    # print curr, end
    
    # res name = 0
    # internal_energy = 1
    # VdW = 4
    # electrostatics = 7
    # polar solvation = 10
    # non-polar solvation = 13
    # total energy = 16
    
    count = 0
    return_dict = {}
    while (curr < end):
        data = lines[curr].split(',')
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict


def parseMMBPSA_backbone_decomp(mmpbsa_file_path, num_res):
    '''

    parse MMPBSA backbone decomposition
    #TODO
    Nick

    Parameters
    ----------
    mmpbsa_file_path
    num_res

    Returns
    -------

    '''
    
    import os
    
    if not os.path.isfile(mmpbsa_file_path):
        return None
    
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Backbone Energy Decomposition" in lines[i]):
            start = i
            break
    
    if (start == -1):
        return None
    
    curr = i + 3
    end = i + 3 + (num_res)
    print(curr, end)
    
    # res name = 0
    # internal_energy = 1
    # VdW = 4
    # electrostatics = 7
    # polar solvation = 10
    # non-polar solvation = 13
    # total energy = 16
    
    count = 0
    return_dict = {}
    while (curr < end):
        data = lines[curr].split(',')
        print((data[0], data[1], data[4], data[7], data[10], data[13], data[16]))
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict


def parseMMBPSA_sidechain_decomp(mmpbsa_file_path, num_res):
    '''

    parse MMPBSA energies from sidechain decomposition

    #TODO
    Nick

    Parameters
    ----------
    mmpbsa_file_path
    num_res

    Returns
    -------

    '''
    import os
    
    if not os.path.isfile(mmpbsa_file_path):
        return None
    
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Sidechain Energy Decomposition" in lines[i]):
            start = i
            break
    
    if (start == -1):
        return None
    
    curr = i + 3
    end = i + 3 + (num_res)
    print(curr, end)
    
    # res name = 0
    # internal_energy = 1
    # VdW = 4
    # electrostatics = 7
    # polar solvation = 10
    # non-polar solvation = 13
    # total energy = 16
    
    count = 0
    return_dict = {}
    while (curr < end):
        data = lines[curr].split(',')
        print((data[0], data[1], data[4], data[7], data[10], data[13], data[16]))
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict


def parseMMPBSA_pairwise_total(mmpbsa_file_path, num_res):
    '''

    parse MMPBSA pairwise decomposition for all residues
    TODO

    Parameters
    ----------
    mmpbsa_file_path
    num_res

    Returns
    -------

    '''
    
    import os
    
    if not os.path.isfile(mmpbsa_file_path):
        return None, None
    
    pair_decomp = np.zeros((num_res, num_res), dtype=np.float)
    
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Total Energy Decomposition" in lines[i]):
            start = i
            break
    
    if (start == -1):
        return None, None
    
    curr = i + 3
    end = i + 3 + (num_res * num_res)
    
    # res name = 0
    # internal_energy = 1
    # VdW = 4
    # electrostatics = 7
    # polar solvation = 10
    # non-polar solvation = 13
    # total energy = 16
    
    count = 0
    return_dict = {}
    while (curr < end):
        data = lines[curr].split(',')
        # print (data[0], data[1], data[4], data[7], data[10], data[13], data[16])

        return_dict[count] = (float(data[2]), float(data[5]), float(data[8]), float(data[11]), float(data[14]), float(data[17]))
        pair_decomp[np.int(old_div(count, num_res))][count % num_res] = float(data[17])

        curr += 1
        count += 1
    # print (pair_decomp)
    return pair_decomp, return_dict


def parseAmberEnergy(amber_file):
    '''

    parses an amber energy file (shown below)

    .. code-block:: none

        #        FINAL RESULTS
        #
        #
        #
        #    NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
        #     100      -2.0167E+03     4.7684E-01     4.3950E+00     C         853
        #
        #  BOND    =       30.0701  ANGLE   =       77.6956  DIHED      =      484.1676
        #  VDWAALS =     -404.1386  EEL     =    -4577.1353  EGB        =     -877.2771
        #  1-4 VDW =      169.7472  1-4 EEL =     3080.1239  RESTRAINT  =        0.0000
        #



    Parameters
    ----------
    amber_file: str
        path to amber.out file.

    Returns
    -------
    final_energy:int
        999 if file not found and exits programm
        finaly energy from amber file
    '''


    if (not os.path.exists(amber_file)):
        print("Amber output file not found")
        return "Error"
    
    amberout = open(amber_file, 'r')
    lines = amberout.readlines()
    amberout.close()
    for i, line in enumerate(lines):
        if ("FINAL RESULTS" in line):
            data = lines[i + 5].split()
    if 'data' in locals():
        if ('NaN' in data[1]):
            return 999
        else:
            return float(data[1])
    else:
        print("Not able to find Amber energy in output file") 
        sys.exit()
    
    return 999


def parseGBEnergy(amber_file):
    """
    parses GB energy output file based on the ff: line containing the total energy of the system from the GB calculation

     .. code-block:: none

            File: _MMPBSA_complex_gb.mdout.X    (x=[0,no_frames])

            ff:     0  -1899.89    788.43   -214.15  -1413.76      0.00  -1060.41  4.75e+00
            BOND    =       35.2852  ANGLE   =      156.4010  DIHED      =      596.7431
            VDWAALS =     -378.8097  EEL     =    -4714.2131  EGB        =    -1060.4105
            1-4 VDW =      164.6609  1-4 EEL =     3300.4508  RESTRAINT  =        0.0000
            ESURF   =        0.0000

    Parameters
    ----------
    amber_file: str
        /path/to/MMPBSAoutput file

    Returns
    -------
    energy: float
        energy value extracted from MMPBSA file or 999 if not found

    """


    if (not os.path.exists(amber_file)):
        print("Amber output file not found")
        return 999

    amberout = open(amber_file, 'r')
    lines = amberout.readlines()

    for i, line in enumerate(lines):
        if ("ff:" in line):
            data = lines[i].split()

    if 'data' in locals():
        if ('NaN' in data[2]):
            return 999
        else:
            return float(data[2])
    else:
        print("Not able to find Amber energy in output file")
        # sys.exit() # commented because sometimes a segfault occurs with pmemd and then this individual is only skipped based on the exit code
    return 999
    
if __name__ == "__main__":
    '''
     Standalone run with first sys argument as the path to the pdb file and test minimization with amber.
    '''
    import sys
    runtleap(mol_file=sys.argv[1])
    runAmberMPI()
     
    fout = open('min_struct.pdb', 'w')
    ferr = open('min_struct.log', 'w')
    try:
        proc = subprocess.Popen(["ambpdb", "-p", "mol.prmtop", "-c", "amber.rst"], stdout=fout, stderr=ferr)
        proc.wait()
        fout.close()
        ferr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("convert rst to pdb failed, returned code %d (check '" + "min_struct.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to convert rst to pdb: %s" % (str(e)))
    runtleap(mol_file="min_struct.pdb")
    
    runMMPBSA()
    
