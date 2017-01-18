'''
 
RotamerCreator.py

@author: Nicholas Browning 
'''
 

#!/software/anaconda-2.0.1/bin/python

import collections
import re
import os
import subprocess
import shutil
import numpy as np

'''
Created on Feb 2, 2016

@author: Nicholas Browning
'''

fail = 999999999

rotamerListAllowed = []
rotamers = []

rotamerEnergyList = collections.OrderedDict()
rotamerPointerList = collections.OrderedDict()

psiphi_data = [[26, 14, 15, 16, 180], [26, 14, 13, 28, 180], [23, 3, 4, 5, 180], [23, 3, 2, 1, 180], [24, 10, 11, 12, 180], [24, 10, 9, 31, 180], [34, 33, 1, 32, 180], [34, 33, 36, 39, 180]]
psiphi_nopro_data = [[26, 14, 15, 16], [26, 14, 13, 28], [24, 10, 11, 12], [24, 10, 9, 31], [34, 33, 1, 32], [34, 33, 36, 39]]


MOL_BUILD_MODE_SUBSTITUTION_ONLY = 0
MOL_BUILD_MODE_DIHEDRAL_ONLY = 1
MOL_BUILD_MODE_SUBST_DIHEDRAL = 2

def buildMolBuildInput(input_structure_path=None, output_path=None, mode=-1, molecule_data=None, frag_data=None, subst_options=None, dihedral_data=None, frag_library="/lcbcdata/nbrownin/Rotamers_mol2", output_structure_path="solution.xyz"):
    output_file = open(output_path, 'w')
    
    output_file.write("%LOAD_MOL\n")
    output_file.write(input_structure_path + "\n\n")
    
    if (mode == MOL_BUILD_MODE_SUBSTITUTION_ONLY or mode == MOL_BUILD_MODE_SUBST_DIHEDRAL):
        output_file.write("%FRAGMENT_LIBRARY\n")
        output_file.write(frag_library + "\n\n")
        
    output_file.write("%OUT_FILE\n")
    output_file.write(output_structure_path + "\n\n")
    
    if (mode == MOL_BUILD_MODE_SUBSTITUTION_ONLY or mode == MOL_BUILD_MODE_SUBST_DIHEDRAL):
        output_file.write("%ACTION\n")
        for i in xrange (0, len(molecule_data)):
            output_file.write("SUBST " + " ".join(str(v) for v in molecule_data[i]) + " " + " ".join(str(v) for v in frag_data[i]) + " " + " ".join(str(v) for v in subst_options[i]) + "\n")
        output_file.write("\n")
    
    if (mode == MOL_BUILD_MODE_DIHEDRAL_ONLY or mode == MOL_BUILD_MODE_SUBST_DIHEDRAL):
        y = len(dihedral_data[0])
        print y
        if (y == 3):
            output_file.write("%SET_PSIPHI\n")
        else:
            output_file.write("%SET_DIHEDRALS\n")
        for i in xrange (0, len(dihedral_data)):
            output_file.write(" ".join(str(v) for v in dihedral_data[i]) + "\n")
    output_file.close()
            
    
def runMolBuild(inputFile="solution.molbuild"):
    output = subprocess.Popen(["./MoleculeBuilder", inputFile], stdout=subprocess.PIPE)
    return output.communicate()

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

def createdata(solution, data, m):
    moleculedata = None
    fragdata = None
    dihedraldata = None
    
    if (m == MOL_BUILD_MODE_SUBST_DIHEDRAL):
        pass
    elif (m == MOL_BUILD_MODE_SUBSTITUTION_ONLY):
        pass
    elif (m == MOL_BUILD_MODE_DIHEDRAL_ONLY):
        dihedraldata = np.zeros((len(solution), 5))
        for i in xrange (len(data)):
            for j in xrange (len(data[i])):
                dihedraldata[i][j] = data[i][j]
            dihedraldata[i][len(dihedraldata[i]) - 1] = solution[i]
    return moleculedata, fragdata, dihedraldata


GAUSSIAN09 = 0

def getBackboneDihedrals(dihedral_data=None, structure_in=None, molbuild_in='print_dihedrals.molbuild'):
    fp = open(molbuild_in, 'w')
    
    fp.write('%LOAD_MOL\n')
    fp.write(structure_in + '\n')
    fp.write('%PRINT_DIHEDRALS\n')
    
    for i in dihedral_data:
        fp.write(" ".join(str(v) for v in i) + '\n')
    fp.close()
    
    output, bla = runMolBuild(inputFile=molbuild_in)
    return np.fromstring(output, sep='\n')

            
# print getBackboneDihedrals(dihedral_data=psiphi_data, structure_in='TOP5_CIS.pdb', molbuild_in='test.mb')


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
            
        
# extractOptimisedGeometry(logfile='TOP5_CIS.log')
# print getBackboneDihedrals(dihedral_data=psiphi_data, structure_in='opt_solution.xyz', molbuild_in='test.mb') 

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
    

# fragdata = [[3 , 5, "M01.mol2"], [3, 5, "M01.mol2"]]
# moldata = [[25, 27], [63 , 65]]
# subst_options = [[1 , 1 , 0 , 1], [1 , 1 , 0 , 1]]
# dihedraldata = [[10, 7, 5, 1, 0], [105, 104, 106, 108, 0]]

# buildMolBuildInput(input_structure_path="test.xyz", output_path="test.molbuild", mode=2, molecule_data=moldata, frag_data=fragdata, subst_options=subst_options, dihedral_data=dihedraldata)
# output = subprocess.Popen(["../MoleculeBuilder/Debug/MoleculeBuilder test.molbuild"], shell=True, stdout=subprocess.PIPE)
# print output.communicate()

# [1, 2, 3, 6 , 0], [5, 4, 3, 6, 0] - PRO

# dihedraldata = [[39, 36, 33, 34, 90], [32, 1, 33, 34, 90], [31, 9, 10, 24, 90], [12, 11, 10, 24, 90], [28, 13, 14, 27, 90], [16, 15, 14, 26, 90]]

# psi = [[34, 33, 36, 39, 0], [24, 10, 9, 31, -90], [27, 14, 13, 28, 90]]
# buildMolBuildInput(input_structure_path="TOP5_CIS.pdb", output_path="test_dhdral.molbuild", mode=MOL_BUILD_MODE_DIHEDRAL_ONLY, dihedral_data=psiphi_nopro)
# output = subprocess.Popen(["../MoleculeBuilder/Debug/MoleculeBuilder test_dhdral.molbuild"], shell=True, stdout=subprocess.PIPE)

# for l in output.communicate():
#    print l


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

def initRotamerLists(dielectric=10, rotamerLibraryPath="/lcbcdata/nbrownin/Rotamers_mol2/"):
    _digits = re.compile('\d')
    
    def contains_digits(d):
        return bool(_digits.search(d))

    rotamerPointerList["ALA"] = 1
    rotamerPointerList["A01"] = 1
        # argenine(ARG)
    rotamerPointerList["ARG"] = 2
    rotamerPointerList["R01"] = 2
    rotamerPointerList["R02"] = 2
    rotamerPointerList["R03"] = 2
    rotamerPointerList["R04"] = 2
    rotamerPointerList["R05"] = 2
    rotamerPointerList["R06"] = 2
    rotamerPointerList["R07"] = 2
    rotamerPointerList["R08"] = 2
    rotamerPointerList["R09"] = 2
    rotamerPointerList["R10"] = 2
    rotamerPointerList["R11"] = 2
    rotamerPointerList["R12"] = 2
    rotamerPointerList["R13"] = 2
    rotamerPointerList["R14"] = 2
    rotamerPointerList["R15"] = 2
    rotamerPointerList["R16"] = 2
    rotamerPointerList["R17"] = 2
    rotamerPointerList["R18"] = 2
    rotamerPointerList["R19"] = 2
    rotamerPointerList["R20"] = 2
    rotamerPointerList["R21"] = 2
    rotamerPointerList["R22"] = 2
    rotamerPointerList["R23"] = 2
    rotamerPointerList["R24"] = 2
    rotamerPointerList["R25"] = 2
    rotamerPointerList["R26"] = 2
    rotamerPointerList["R27"] = 2
    rotamerPointerList["R28"] = 2
        # hydrogenated aspartic acid (ASH)
    rotamerPointerList["ASH"] = 3
    rotamerPointerList["D06"] = 3
    rotamerPointerList["D07"] = 3
    rotamerPointerList["D08"] = 3
    rotamerPointerList["D09"] = 3
    rotamerPointerList["D10"] = 3
        # aspartic acid (ASP)
    rotamerPointerList["ASP"] = 4
    rotamerPointerList["D01"] = 4
    rotamerPointerList["D02"] = 4
    rotamerPointerList["D03"] = 4
    rotamerPointerList["D04"] = 4
    rotamerPointerList["D05"] = 4
        # asparagine (ASN)
    rotamerPointerList["ASN"] = 5
    rotamerPointerList["N01"] = 5
    rotamerPointerList["N02"] = 5
    rotamerPointerList["N03"] = 5
    rotamerPointerList["N04"] = 5
    rotamerPointerList["N05"] = 5
    rotamerPointerList["N06"] = 5
    rotamerPointerList["N07"] = 5
        # Cysteine (CYS)
    rotamerPointerList["CYS"] = 6
    rotamerPointerList["C01"] = 6
    rotamerPointerList["C02"] = 6
    rotamerPointerList["C03"] = 6
        # hydrogenated glutamic acid (GLH)
    rotamerPointerList["GLH"] = 7
    rotamerPointerList["E08"] = 7
    rotamerPointerList["E09"] = 7
    rotamerPointerList["E10"] = 7
    rotamerPointerList["E11"] = 7
    rotamerPointerList["E12"] = 7
    rotamerPointerList["E13"] = 7
    rotamerPointerList["E14"] = 7
        # glutamic acid (GLU)
    rotamerPointerList["GLU"] = 8
    rotamerPointerList["E01"] = 8
    rotamerPointerList["E02"] = 8
    rotamerPointerList["E03"] = 8
    rotamerPointerList["E04"] = 8
    rotamerPointerList["E05"] = 8
    rotamerPointerList["E06"] = 8
    rotamerPointerList["E07"] = 8
        # glutamine (GLN)
    rotamerPointerList["GLN"] = 9
    rotamerPointerList["Q01"] = 9
    rotamerPointerList["Q02"] = 9
    rotamerPointerList["Q03"] = 9
    rotamerPointerList["Q04"] = 9
    rotamerPointerList["Q05"] = 9
        # glycine (GLY)
    rotamerPointerList["GLY"] = 10
    rotamerPointerList["G01"] = 10
        # HID
    rotamerPointerList["HID"] = 11
    rotamerPointerList["HD1"] = 11
    rotamerPointerList["HD2"] = 11
    rotamerPointerList["HD3"] = 11
    rotamerPointerList["HD4"] = 11
    rotamerPointerList["HD5"] = 11
    rotamerPointerList["HD6"] = 11
    rotamerPointerList["HD7"] = 11
    rotamerPointerList["HD8"] = 11
        # HIE
    rotamerPointerList["HIE"] = 12
    rotamerPointerList["HE1"] = 12
    rotamerPointerList["HE2"] = 12
    rotamerPointerList["HE3"] = 12
    rotamerPointerList["HE4"] = 12
    rotamerPointerList["HE5"] = 12
    rotamerPointerList["HE6"] = 12
    rotamerPointerList["HE7"] = 12
    rotamerPointerList["HE8"] = 12
        # HIP
    rotamerPointerList["HIP"] = 13
    rotamerPointerList["HP1"] = 13
    rotamerPointerList["HP2"] = 13
    rotamerPointerList["HP3"] = 13
    rotamerPointerList["HP4"] = 13
    rotamerPointerList["HP5"] = 13
    rotamerPointerList["HP6"] = 13
    rotamerPointerList["HP7"] = 13
    rotamerPointerList["HP8"] = 13
        # isoleucine (ILE)
    rotamerPointerList["ILE"] = 14
    rotamerPointerList["I01"] = 14
    rotamerPointerList["I02"] = 14
    rotamerPointerList["I03"] = 14
    rotamerPointerList["I04"] = 14
    rotamerPointerList["I05"] = 14
        # Leucine (LEU)
    rotamerPointerList["LEU"] = 15
    rotamerPointerList["L01"] = 15
    rotamerPointerList["L02"] = 15
    rotamerPointerList["L03"] = 15
    rotamerPointerList["L04"] = 15
        # lysine (LYN)
    rotamerPointerList["LYN"] = 16
    rotamerPointerList["K20"] = 16
    rotamerPointerList["K21"] = 16
    rotamerPointerList["K22"] = 16
    rotamerPointerList["K23"] = 16
    rotamerPointerList["K24"] = 16
    rotamerPointerList["K25"] = 16
    rotamerPointerList["K26"] = 16
    rotamerPointerList["K27"] = 16
    rotamerPointerList["K28"] = 16
    rotamerPointerList["K29"] = 16
    rotamerPointerList["K30"] = 16
    rotamerPointerList["K31"] = 16
    rotamerPointerList["K32"] = 16
    rotamerPointerList["K33"] = 16
    rotamerPointerList["K34"] = 16
    rotamerPointerList["K35"] = 16
    rotamerPointerList["K36"] = 16
    rotamerPointerList["K37"] = 16
    rotamerPointerList["K38"] = 16
        # lysine (LYS)
    rotamerPointerList["LYS"] = 17
    rotamerPointerList["K01"] = 17
    rotamerPointerList["K02"] = 17
    rotamerPointerList["K03"] = 17
    rotamerPointerList["K04"] = 17
    rotamerPointerList["K05"] = 17
    rotamerPointerList["K06"] = 17
    rotamerPointerList["K07"] = 17
    rotamerPointerList["K08"] = 17
    rotamerPointerList["K09"] = 17
    rotamerPointerList["K10"] = 17
    rotamerPointerList["K11"] = 17
    rotamerPointerList["K12"] = 17
    rotamerPointerList["K13"] = 17
    rotamerPointerList["K14"] = 17
    rotamerPointerList["K15"] = 17
    rotamerPointerList["K16"] = 17
    rotamerPointerList["K17"] = 17
    rotamerPointerList["K18"] = 17
    rotamerPointerList["K19"] = 17
        # methionine (MET)
    rotamerPointerList["MET"] = 18
    rotamerPointerList["M01"] = 18
    rotamerPointerList["M02"] = 18
    rotamerPointerList["M03"] = 18
    rotamerPointerList["M04"] = 18
    rotamerPointerList["M05"] = 18
    rotamerPointerList["M06"] = 18
    rotamerPointerList["M07"] = 18
    rotamerPointerList["M08"] = 18
    rotamerPointerList["M09"] = 18
    rotamerPointerList["M10"] = 18
    rotamerPointerList["M11"] = 18
    rotamerPointerList["M12"] = 18
    rotamerPointerList["M13"] = 18
        # (PHE)
    rotamerPointerList["PHE"] = 19
    rotamerPointerList["P01"] = 19
    rotamerPointerList["P02"] = 19
    rotamerPointerList["P03"] = 19
    rotamerPointerList["P04"] = 19
        # serine (SER)
    rotamerPointerList["SER"] = 20
    rotamerPointerList["S01"] = 20
    rotamerPointerList["S02"] = 20
    rotamerPointerList["S03"] = 20
        # threonine (THR)
    rotamerPointerList["THR"] = 21
    rotamerPointerList["T01"] = 21
    rotamerPointerList["T02"] = 21
    rotamerPointerList["T03"] = 21
        # tryptophan (TRP)
    rotamerPointerList["TRP"] = 22
    rotamerPointerList["W01"] = 22
    rotamerPointerList["W02"] = 22
    rotamerPointerList["W03"] = 22
    rotamerPointerList["W04"] = 22
    rotamerPointerList["W05"] = 22
    rotamerPointerList["W06"] = 22
    rotamerPointerList["W07"] = 22
        # tyrosine (TYR)
    rotamerPointerList["TYR"] = 23
    rotamerPointerList["Y01"] = 23
    rotamerPointerList["Y02"] = 23
    rotamerPointerList["Y03"] = 23
    rotamerPointerList["Y04"] = 23
        # valine (VAL)
    rotamerPointerList["VAL"] = 24
    rotamerPointerList["V01"] = 24
    rotamerPointerList["V02"] = 24
    rotamerPointerList["V03"] = 24
    # PRO
    rotamerPointerList["PRO"] = 25
    
    
    if dielectric == 10:
        rotamerEnergyList[1] = 13.03
        rotamerEnergyList[2] = -170.77
        rotamerEnergyList[3] = -43.74
        rotamerEnergyList[4] = -52.27
        rotamerEnergyList[5] = -72.72
        rotamerEnergyList[6] = 13.55
        rotamerEnergyList[7] = -36.86
        rotamerEnergyList[8] = -47.00
        rotamerEnergyList[9] = -56.93
        rotamerEnergyList[10] = 7.11
        rotamerEnergyList[11] = 12.43
        rotamerEnergyList[12] = 8.74
        rotamerEnergyList[13] = 22.66
        rotamerEnergyList[14] = 10.74
        rotamerEnergyList[15] = -10.04
        rotamerEnergyList[16] = -4.36
        rotamerEnergyList[17] = -9.60
        rotamerEnergyList[18] = 8.36
        rotamerEnergyList[19] = 11.05
        rotamerEnergyList[20] = -0.71
        rotamerEnergyList[21] = -19.37
        rotamerEnergyList[22] = 14.43
        rotamerEnergyList[23] = -13.28
        rotamerEnergyList[24] = -7.58
        rotamerEnergyList[25] = 0
    elif dielectric == 50:
        rotamerEnergyList[1] = 12.14
        rotamerEnergyList[2] = -176.22
        rotamerEnergyList[3] = -45.64
        rotamerEnergyList[4] = -58.92
        rotamerEnergyList[5] = -74.29
        rotamerEnergyList[6] = 12.58
        rotamerEnergyList[7] = -38.93
        rotamerEnergyList[8] = -53.61
        rotamerEnergyList[9] = -58.79
        rotamerEnergyList[10] = 6.15
        rotamerEnergyList[11] = 10.83
        rotamerEnergyList[12] = 7.22
        rotamerEnergyList[13] = 17.17
        rotamerEnergyList[14] = 9.93
        rotamerEnergyList[15] = -10.88
        rotamerEnergyList[16] = -5.73
        rotamerEnergyList[17] = -15.72
        rotamerEnergyList[18] = 7.47
        rotamerEnergyList[19] = 10.04
        rotamerEnergyList[20] = -2.32
        rotamerEnergyList[21] = -20.45
        rotamerEnergyList[22] = 13.07
        rotamerEnergyList[23] = -14.89
        rotamerEnergyList[24] = -8.41
        rotamerEnergyList[25] = 0
    elif dielectric == 80:
        rotamerEnergyList[1] = 12.06
        rotamerEnergyList[2] = -176.73
        rotamerEnergyList[3] = -45.82
        rotamerEnergyList[4] = -59.55
        rotamerEnergyList[5] = -74.44
        rotamerEnergyList[6] = 12.49
        rotamerEnergyList[7] = -39.13
        rotamerEnergyList[8] = -54.30
        rotamerEnergyList[9] = -58.97
        rotamerEnergyList[10] = 6.06
        rotamerEnergyList[11] = 10.68
        rotamerEnergyList[12] = 7.08
        rotamerEnergyList[13] = 16.62
        rotamerEnergyList[14] = 9.85
        rotamerEnergyList[15] = -10.96
        rotamerEnergyList[16] = -5.85
        rotamerEnergyList[17] = -16.29
        rotamerEnergyList[18] = 7.29
        rotamerEnergyList[19] = 9.95
        rotamerEnergyList[20] = -2.47
        rotamerEnergyList[21] = -20.59
        rotamerEnergyList[22] = 12.94
        rotamerEnergyList[23] = -15.05
        rotamerEnergyList[24] = -8.49
        rotamerEnergyList[25] = 0
    
    rotamersAll = []
    for key, value in rotamerPointerList.items():
        rotamersAll.append(key)
        
    
    for x in xrange (0, len(rotamersAll)):
        if contains_digits(rotamersAll[x]):
            rotamers.append(rotamersAll[x])

    rotamerListFullLibrary = []

    for root, dirs, files in os.walk(rotamerLibraryPath):
        for file in files:
            if ".mol2" in file:
                rotamerListFullLibrary.append(root + "/" + file)
    
    for i in range(len(rotamers)):
        for j in range(len(rotamerListFullLibrary)):
            if rotamers[i] in rotamerListFullLibrary[j]:
                rotamerListAllowed.append(rotamerListFullLibrary[j])
                break
            
    if len(rotamers) != len(rotamerListAllowed):
        print -1, "Problem with rotamer library size!", len (rotamers), len (rotamerListAllowed)
 
  
# createPDB("250ns.pdb", "test.pdb", [5, 36], [0, 0])
# print evaluateRotamerLEAPGeneration("250ns.pdb", originalEnergy=-2016.7, moleculeSubstitutionIndexes=[22, 79, 98], designVector=[50, 45, 51])
