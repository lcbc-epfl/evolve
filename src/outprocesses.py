'''
Created on Jan 10, 2018

@author: Nicholas Browning & Justin Villard

All functions used to run softwares out of the GA code
'''


import os
import subprocess

def runtleap(work_dir="", mol_file="mol.pdb", tleaptemp="leap_template.in", tleapin="leap.in", inpcrd_out="mol.inpcrd", prmtop_out="mol.prmtop"):
    # TODO
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
    
    i = open(work_dir + tleapin, 'w')
    for line in g:
        i.write(line)
    i.close()

    # proc = subprocess.Popen([tleap_path+"/bin/tleap", "-f", work_dir+tleapin], stdout=subprocess.PIPE)
    # out, err = proc.communicate()

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

    
def runAmberMPI(work_dir="", np=2, amberin="amber_minimization.in", prmtop="mol.prmtop", inpcrd="mol.inpcrd", amberout="amber.out", restart="amber.rst"):

    ''' Actually already made by amber with flag -O
    if (os.path.exists(work_dir + amberout)):
        os.remove(work_dir + amberout)
    if (os.path.exists(work_dir + restart)):
        os.remove(work_dir + restart)
    '''    
    # proc = subprocess.Popen(["mpirun", "-np", str(np), "sander.MPI", "-O", "-i", amberin, "-o", work_dir+amberout, "-c", work_dir+inpcrd, "-p", work_dir+prmtop, "-r", work_dir+restart], stdout=subprocess.PIPE)
    # can use pmemd if not in vaccum, more efficient parallelization but apparently doesn't work in vaccum
    # out, err = proc.communicate()

    foutanderr = open(work_dir + '/amber.log', 'w')
    try:
        proc = subprocess.Popen(["mpirun", "-np", str(np), "sander.MPI", "-O", "-i", amberin, "-o", work_dir + amberout, "-c", work_dir + inpcrd, "-p", work_dir + prmtop, "-r", work_dir + restart], stdout=foutanderr, stderr=subprocess.STDOUT)
        proc.wait()
        foutanderr.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except subprocess.CalledProcessError as e:
        sys.exit("sander failed, returned code %d (check '" + work_dir + "/amber.log')" % (e.returncode))
    except OSError as e:
        sys.exit("failed to execute sander.MPI: %s" % (str(e)))
    

def runMMPBSA(workdir="", mmpbsa_in="mmbpsa.in", prmtop="mol.prmtop", inpcrd="mol.inpcrd", final_results="mmpbsa.dat", decomp_results="mmpbsa_decomp.dat"):
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

def parseMMBPSA_total_decomp(mmpbsa_file_path, num_res):
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Total Energy Decomposition" in lines[i]):
            start = i
            break
    
    
    curr = i + 3
    end = i + 3 + (num_res)
    #print curr, end
    
    
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
        #print (data[0], data[1], data[4], data[7], data[10], data[13], data[16])
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict

def parseMMBPSA_backbone_decomp(mmpbsa_file_path, num_res):
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Backbone Energy Decomposition" in lines[i]):
            start = i
            break
    
    
    curr = i + 3
    end = i + 3 + (num_res)
    print curr, end
    
    
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
        print (data[0], data[1], data[4], data[7], data[10], data[13], data[16])
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict

def parseMMBPSA_sidechain_decomp(mmpbsa_file_path, num_res):
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Sidechain Energy Decomposition" in lines[i]):
            start = i
            break
    
    
    curr = i + 3
    end = i + 3 + (num_res)
    print curr, end
    
    
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
        print (data[0], data[1], data[4], data[7], data[10], data[13], data[16])
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict

def parseAmberEnergy(amber_file):
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
    # --------------------------------------------------------------------------------
    #    5.  TIMINGS
    # --------------------------------------------------------------------------------
    #
    if (not os.path.exists(amber_file)):
        print("Amber output file not found")
        return "Error"
    
    amberout = open(amber_file, 'r')
    lines = amberout.readlines()

    for i, line in enumerate(lines):
        if ("FINAL RESULTS" in line):
            data = lines[i + 5].split()

    if 'data' in locals():
        return float(data[1])
    else:
        print("Not able to find Amber energy in output file") 
        raise EOFError
