
'''
Created on Feb 2, 2016

@author: Nicholas Browning
'''

import collections
import re
import os
import subprocess
import shutil
import numpy as np

import openbabel

fail = 999999999

def buildg09InputFile(templateFile='template.com', molbuildoutputstructure='solution.xyz', outputFile='solution.com'):
    
    templateInputFile = open(templateFile, 'r')
    
    templateInputFileLines = templateInputFile.readlines()
    
    templateInputFile.close()
    
    outputFile = open(outputFile, 'w')
    
    for v in templateInputFileLines:
        outputFile.write(v)
        
    molbuild_output = open (molbuildoutputstructure, 'r')
    
    molbuild_output_lines = molbuild_output.readlines()
    
    molbuild_output.close()
    
    for i in xrange (2, len(molbuild_output_lines)):
        outputFile.write(molbuild_output_lines[i])
    outputFile.write("\n\n\n")
    outputFile.close()

def extractg09Energy(filein='solution.log'):
    # SCF Done:  E(RAM1) = -0.193407985355     A.U. after    8 cycles
    #       NFock=  7  Conv=0.72D-08     -V/T= 0.9994
    # SCF Done:  E(UHF) =  -1019.20269909     A.U. after   11 cycles
    
    fp = open(filein, 'r')
    
    lines = fp.readlines()
    
    for i in xrange (len(lines) - 1, 0, -1):
        if ('SCF Done:' in lines[i]):
            splitted = lines[i].split()
            return float(splitted[4])
            
    return 9999999

def extractOptimisedGeometry(program=GAUSSIAN09, out_structure='opt_solution.xyz', **kwargs):
    if (program == GAUSSIAN09):
        fchk_file = kwargs.setdefault('fchkfile', None)
        log_file = kwargs.setdefault('logfile', None)
        if (log_file != None):
            output = subprocess.Popen(["babel", "-ig09", log_file, "-oxyz", out_structure], stdout=subprocess.PIPE)
            return output.communicate()
        else:
            output = subprocess.Popen(["babel", "-ifchk", fchk_file, "-oxyz", out_structure], stdout=subprocess.PIPE)
            return output.communicate()
            

def evaluateRotamerMolBuilderGaussian(initialMoleculePath=None, gaussianTemplateCom="template.com", gaussianOutputName="solution", outputMolbuildFile="solution.molbuild", outputStructurePath='solution.xyz', designVector=None, m=None, data=None):

    moleculedata = None
    fragdata = None
    dihedraldata = None
    
    if (m == MOL_BUILD_MODE_SUBST_DIHEDRAL):
        pass
    elif (m == MOL_BUILD_MODE_SUBSTITUTION_ONLY):
        pass
    elif (m == MOL_BUILD_MODE_DIHEDRAL_ONLY):
        dihedraldata = np.zeros((len(designVector), 5))
        for i in xrange (len(data)):
            for j in xrange (len(data[i])):
                dihedraldata[i][j] = data[i][j]
            dihedraldata[i][len(dihedraldata[i]) - 1] = designVector[i]

            
    buildMolBuildInput(input_structure_path=initialMoleculePath, output_path=outputMolbuildFile, output_structure_path=outputStructurePath, mode=m, molecule_data=moleculedata, frag_data=fragdata, dihedral_data=dihedraldata)
    runMolBuild(inputFile=outputMolbuildFile)
    buildg09InputFile(templateFile=gaussianTemplateCom, molbuildoutputstructure=outputStructurePath, outputFile=gaussianOutputName + '.com')
    p = subprocess.Popen(["g09", gaussianOutputName + '.com'], stdout=subprocess.PIPE)
    p.communicate()
    return extractg09Energy(gaussianOutputName + '.log')
    
def writeOBMol(obmol, filename):
    obconv = openbabel.OBConversion()
    
    if (not obconv.setOutFormat(filename.split(".")[1])):
        print "Could not set output file format"
        return False
    
    return obconv.WriteFile(obmol, filename)

def evaluateBackboneFitness_amber(settings, mol, orig_energy):
    if (not writeOBMol(mol, "solution.pdb")):
        return None
    
    runtleap(tleapin=settings.tleap_template)
    runAmberMPI(np=settings.MPI_procs)
    
    finalEnergy = parseAmberEnergy()
    
    if (finalEnergy == fail):
        return fail
    
    

def evaluateRotamerLEAPGeneration(initialMoleculePath, originalEnergy, moleculeSubstitutionIndexes, designVector):
    tempMolName = "tmp.pdb"
    shutil.copy(initialMoleculePath, tempMolName);
    resChanges = createPDB(tempMolName, "solution.pdb", moleculeSubstitutionIndexes, designVector)
    print resChanges
    runtleap()
    runAmberMPI(np=8)
    
    finalEnergy = parseAmberEnergy()
    
    if (finalEnergy == fail):
        return fail
            
    add = 0.0;
    negate = originalEnergy;
    
    addValues = []
    negateValues = []

    for i in range (len(resChanges)):
    
        sub = resChanges[i].split("->")
    
        if sub[0] != sub[1]:
            addValue = rotamerEnergyList[rotamerPointerList[sub[0][:3]]]
            addValues.append(addValue)
            negateValue = rotamerEnergyList[rotamerPointerList[sub[1][:3]]]
            negateValues.append(negateValue)
            add += addValue;
            negate += negateValue;
    
    print add, negate, finalEnergy, addValues, negateValues        
    finalEnergy += (add - negate)
    # print finalEnergy
    # print resChanges, "add: ", add, "negate: ", negate
    return finalEnergy


def runtleap(pdb="solution.pdb", tleapin="tleap.in", tleapout="tleap_ga.in", inpcrd_out="solution.inpcrd", prmtop_out="solution.prmtop"):
    # TODO
    if (os.path.exists(inpcrd_out)):
        os.remove(inpcrd_out)
    if (os.path.exists(prmtop_out)):
        os.remove(prmtop_out)
    if (os.path.exists("leap.log")):
        os.remove("leap.log")

    f = open(tleapin, "r")
    g = f.readlines()
    
    for i, line in enumerate(g):
        if "solution=" in line:
            g[i] = "solution=loadpdb \"" + pdb + "\"\n"
        elif "saveamberparm" in line:
            g[i] = "saveamberparm solution " + prmtop_out + " " + inpcrd_out + "\n"
    
    i = open(tleapout, 'w')
    for line in g:
        i.write(line)
    i.close()
    
    proc = subprocess.Popen(["/software/amber/amber12/bin/tleap", "-f", tleapout], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    
  
def runAmberMPI(np, prmtop="solution.prmtop", inpcrd="solution.inpcrd", amber_out="amber_out.out" , amberin="amber.in", restart="amber_out.rst"):

    if (os.path.exists(amber_out)):
        os.remove(amber_out)
    if (os.path.exists(restart)):
        os.remove(restart)
        
    proc = subprocess.Popen(["mpirun", "-np", str(np), "sander.MPI", "-O", "-i", amberin, "-o", amber_out, "-c", inpcrd, "-p", prmtop, "-r", restart ], stdout=subprocess.PIPE)
    out, err = proc.communicate()


'''Parses ENERGY in FINAL RESULTS section'''
def parseAmberEnergy(amber_out="amber_out.out"):
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
    if (not os.path.exists(amber_out)):
        return fail
    
    amberout = open(amber_out, 'r')
   
    lines = amberout.readlines()
    
    for i, v in enumerate(lines):
        if ("FINAL RESULTS" in v):
            data = lines[i + 5].split()
            return float(data[1])
            break
    return fail


def createPDB(modifyPDBPath, finalPDBPath, moleculeSubstitutionIndexes, designVector):
    resChanges = []
    readFile = open(modifyPDBPath, 'r')
    
    lines = readFile.readlines()
    
    readFile.close()
    
    cutStartIndexes = []
    cutFinishIndexes = []
    startMoleculeIndexes = []
    endMoleculeIndexes = []
    mutateRes = []
    
    for i, v in enumerate (designVector):
        index = int(v)
        startPosition = moleculeSubstitutionIndexes[i]

       
        endMolecule = -1
        startMolecule = -1

        current = startPosition
        CA_lineData = lines[current].split()
        CA_resName = CA_lineData[3].strip()
        resChanges.append(CA_resName + "->" + rotamers[index])
        mutateRes.append(rotamers[index])
        CA_resNum = int(CA_lineData[4])

        # forward scan
        while (True):
            
            if (current > len(lines) - 1):
                endMolecule = current - 1
                break

            line = lines[current]
            lineData = line.split()
            
            if ("ATOM" not in lineData):
                endMolecule = current - 1
                break
        
            resNum = int(lineData[4])
            
            if (resNum != CA_resNum):
                endMolecule = current - 1
                break
            current += 1
            
        current = startPosition
        
        # backward scan
        while (True):
            
            if (current < 0):
                startMolecule = current + 1
                break

            line = lines[current]
            lineData = line.split()
            
            if ("ATOM" not in lineData):
                startMolecule = current + 1
                break
        
            resNum = int(lineData[4])
            
            if (resNum != CA_resNum):
                startMolecule = current + 1
                break
            
            current -= 1
            
        # find C terminus
        
        cterminus = -1
        for j in xrange(startMolecule, endMolecule):
            line = lines[j]
            lineData = line.split()
            atomName = lineData[2].strip()
            
            if (atomName == "C"):
                cterminus = j
                break
        
        if (startMolecule == -1 or endMolecule == -1 or cterminus == -1):
            print "ERROR: startMol,endMol or cterminus ==-1!"
            return 
            
                
        cutFinishIndexes.append(cterminus)    
        cutStartIndexes.append(startPosition)
        startMoleculeIndexes.append(startMolecule)
        endMoleculeIndexes.append(endMolecule)
    
    outFile = open(finalPDBPath, 'w')
    
    for i in xrange(0, len(startMoleculeIndexes)):
        
        currStartMol = startMoleculeIndexes[i]
        currEndMol = endMoleculeIndexes[i]
        curCutStart = cutStartIndexes[i]
        currCutFinish = cutFinishIndexes[i]
        startLineData = lines[currStartMol].split()
    
        currResName = startLineData[3]
        
        if ("PRO" not in currResName):
            for j in xrange (currStartMol, currEndMol + 1):
                line = lines[j]
                if (j >= curCutStart + 2 and j <= currCutFinish - 1):
                    # print lines[j]
                    lines[j] = None
                else:
                    lines[j] = line.replace(line[17:20], mutateRes[i])
        else:
            for j in xrange (currStartMol, currEndMol + 1):
                line = lines[j]
                if (j == currStartMol or j == currEndMol):
                    lines[j] = line.replace(line[17:20], mutateRes[i])
                else:
                    lines[j] = None

      
    for i in xrange(0, len(lines)):
        if (lines[i] is not None):
            outFile.write(lines[i])
            
    outFile.close()
    return resChanges
