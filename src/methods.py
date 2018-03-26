'''
Created on Nov 4, 2015

@author: roethlisbergergroup
'''


import math
import numpy


element2mass = {'Ru': 101.904348, 'Re': 186.955744, 'Rf': 261.10869, 'Ra': 226.025403,
                'Rb': 84.9118, 'Rn': 222.017571, 'Rh': 102.905503, 'Be': 9.012183,
                'Ba': 137.905232, 'Bh': 262.12293, 'Bi': 208.980374, 'Bk': 247.0703,
                'Br': 78.918336, 'H': 1.007825, 'P': 30.973763, 'Os': 191.961467, 'Es': 252.082944,
                'Hg': 201.970617, 'Ge': 73.921179, 'Gd': 157.924019, 'Ga': 68.925581,
                'Pr': 140.907647, 'Pt': 194.964766, 'Pu': 244.064199, 'C': 12.0,
                'Pb': 207.976627, 'Pa': 231.03588, 'Pd': 105.903475, 'Cd': 113.903361,
                'Po': 208.982404, 'Pm': 144.912743, 'Hs': 269.1341, 'Ho': 164.930319,
                'Hf': 179.946546, 'K': 38.963708, 'He': 4.002603, 'Md': 258.09857,
                'Mg': 23.985045, 'Mo': 97.905405, 'Mn': 54.938046, 'O': 15.994915,
                'Mt': 267.138, 'S': 31.972072, 'W': 183.950928, 'Zn': 63.929145, 'Eu': 152.921225,
                'Zr': 89.904708, 'Er': 165.93029, 'Ni': 57.935347, 'No': 259.100931,
                'Na': 22.98977, 'Nb': 92.906378, 'Nd': 141.907719, 'Ne': 19.992439,
                'Np': 237.048168, 'Fr': 223.019733, 'Fe': 55.934939, 'Fm': 257.095099,
                'B': 11.009305, 'F': 18.998403, 'Sr': 87.905625, 'N': 14.003074,
                'Kr': 83.911506, 'Si': 27.976928, 'Sn': 119.902199, 'Sm': 151.919728,
                'V': 50.943963, 'Sc': 44.955914, 'Sb': 120.903824, 'Sg': 263.11822,
                'Se': 79.916521, 'Co': 58.933198, 'Cm': 247.070347, 'Cl': 34.968853,
                'Ca': 39.962591, 'Cf': 251.07958, 'Ce': 139.905433, 'Xe': 131.904148,
                'Lu': 174.94077, 'Cs': 132.905429, 'Cr': 51.94051, 'Cu': 62.929599,
                'La': 138.906347, 'Li': 7.016005, 'Tl': 204.974401, 'Tm': 168.934212,
                'Lr': 260.10532, 'Th': 232.038051, 'Ti': 47.947947, 'Te': 129.906229,
                'Tb': 158.925342, 'Tc': 97.907216, 'Ta': 180.947992, 'Yb': 173.938859,
                'Db': 262.11376, 'Dy': 163.929171, 'I': 126.904477, 'U': 238.050785,
                'Y': 88.905856, 'Ac': 227.02775, 'Ag': 106.905095, 'Ir': 192.962917,
                'Am': 243.061373, 'Al': 26.981541, 'As': 74.921596, 'Ar': 39.962383,
                'Au': 196.966543, 'At': 209.987126, 'In': 114.903875}

element2charge = {'H':1, 'C': 6, 'N':7, 'O':8, 'F':9, 'Pb':82, 'I':53, 'Cs':55 }
bohr2angstrom = 0.5291772083
angstrom2bohr = 1.88972613386246


#===============================================================================
# 5
# 000001
# C    -0.0126981359     1.0858041578     0.0080009958    -0.535689
# H     0.002150416    -0.0060313176     0.0019761204     0.133921
# H     1.0117308433     1.4637511618     0.0002765748     0.133922
# H    -0.540815069     1.4475266138    -0.8766437152     0.133923
# H    -0.5238136345     1.4379326443     0.9063972942     0.133923
#===============================================================================

def loadXYZ(filename, ang2bohr=True):
    with open(filename, 'r') as f:
        lines = f.readlines()
        numAtoms = int(lines[0])
        positions = numpy.zeros((numAtoms, 3), dtype=numpy.float64)
        elems = [None] * numAtoms
        comment = lines[1]
        for x in xrange (2, 2 + numAtoms):
            line_split = lines[x].rsplit()
            elems[x - 2] = line_split[0]
            positions[x - 2][0] = float(line_split[1]) 
            positions[x - 2][1] = float(line_split[2]) 
            positions[x - 2][2] = float(line_split[3])
            if (ang2bohr):
                positions[x - 2][0] *= angstrom2bohr
                positions[x - 2][1] *= angstrom2bohr
                positions[x - 2][2] *= angstrom2bohr
                
    return elems, positions, comment


def writeXYZ(fileName, elems, positions, bohr2ang):
    with open (fileName, 'w') as f:
        f.write(str(len(elems)) + "\n")
        f.write(" " + "\n")
        for x in xrange (0, len(elems)):
            f.write(elems[x] + " " + " ".join(str(v) for v in positions[x]) + "\n")
    f.close()


def generateRandomTrainingSet(size, trainingDat):
    trainingSet = numpy.empty((size), dtype=numpy.int32)
    allowedTrainingMols = numpy.loadtxt(trainingDat, dtype=numpy.int32, delimiter='\n')
    trainingSet.fill(-1)
    for i in xrange(size):
        val = allowedTrainingMols[numpy.random.randint(0, len(allowedTrainingMols))]
        while (val in trainingSet):
            val = allowedTrainingMols[numpy.random.randint(0, len(allowedTrainingMols))]
        trainingSet[i] = val
    return trainingSet


def createCoulombMatrixNaive(elems, positions, cmsize):
    numAtoms = len(elems)
    cm = numpy.zeros((cmsize, cmsize), dtype=numpy.float64)
    for i in xrange(0, numAtoms):
        for j in xrange (0, numAtoms):
            if (i == j):
                cm[i][j] = 0.5 * math.pow(element2charge[elems[i]], 2.4)
            else:
                dR = numpy.subtract(positions[i], positions[j]) 
                #
                # print (i, j, dR, numpy.dot(dR, dR), numpy.sqrt(numpy.dot(dR, dR)),   numpy.multiply(element2charge[elems[i]], element2charge[elems[j]]), numpy.divide(numpy.multiply(element2charge[elems[i]], element2charge[elems[j]]), numpy.sqrt(numpy.dot(dR, dR))))
                cm[i][j] = (element2charge[elems[i]] * element2charge[elems[j]]) / numpy.sqrt(numpy.dot(dR, dR))
    return cm

def createCoulombMatrix(elems, positions, cmsize):
    numAtoms = len(elems)
    cm = numpy.zeros((cmsize, cmsize), dtype=np.float32)
    for i in xrange(0, numAtoms):
        cm[i][i] = 0.5 * math.pow(element2charge[elems[i]], 2.4)
        for j in xrange (i + 1, numAtoms):
            dR = positions[i] - positions[j]
            cm[i][j] = (element2charge[elems[i]] * element2charge[elems[j]]) / numpy.sqrt(numpy.dot(dR, dR))
            cm[j][i] = cm[i][j]
    return cm

def createRIJMatrix(elems, positions, cmsize):
    numAtoms = len(elems)
    cm = numpy.zeros((cmsize, cmsize), dtype=numpy.float32)
    for i in xrange(0, numAtoms):
        cm[i][i] = 0
        for j in xrange (i + 1, numAtoms):
            dR = positions[i] - positions[j]
            cm[i][j] = 1 / numpy.linalg.norm(dR)
            cm[j][i] = cm[i][j]
    return cm

def CoM(elems, xyzs):
    v = numpy.zeros(3)
    M = 0
    for i in xrange (0, len(elems)):
        v = v + numpy.multiply(element2mass[elems[i]], xyzs[i])
        M = M + element2mass[elems[i]]
    v = v / M
    return v


def eucdistance (cm1, cm2):
    indices = numpy.triu_indices(len(cm1))
    cm1_1 = cm1[indices]
    cm2_1 = cm2[indices]
    return numpy.sqrt(numpy.sum(numpy.square(cm1_1 - cm2_1)))

def distance (cm1, cm2):
    indices = numpy.triu_indices(len(cm1))
    cm1_1 = cm1[indices]
    cm2_1 = cm2[indices]
    return numpy.sum(numpy.abs(cm1_1 - cm2_1))

def createSortedCoulombMatrix(elems, positions, cmsize):
    numAtoms = len(elems)
    cm = numpy.zeros((cmsize, cmsize))
    for i in xrange(0, numAtoms):
        cm[i][i] = 0.5 * math.pow(element2charge[elems[i]], 2.4)
        for j in xrange (i + 1, numAtoms):
            dR = (positions[i] - positions[j])
            cm[i][j] = (element2charge[elems[i]] * element2charge[elems[j]]) / numpy.sqrt(numpy.dot(dR, dR))
            cm[j][i] = cm[i][j]
            
    norms = numpy.zeros((cmsize))
    indexes = numpy.zeros((cmsize))
    for i in xrange (0, cmsize):
        norms[i] = numpy.sqrt(numpy.dot(cm[i], cm[i]))
        indexes[i] = i
        
      
    cm2 = numpy.zeros((cmsize, cmsize))
    for i in xrange (0, cmsize):
        cm2[i][i] = cm[i][i]
        for j in xrange (i + 1, cmsize):
            cm2[i][j] = cm[i][j]
            cm2[j][i] = cm2[i][j]
    
    for i in xrange (0, cmsize):
        for j in xrange (i + 1, cmsize):
            if (norms[i] < norms [j]):
                i1 = indexes[i]
                indexes[i] = indexes[j]
                indexes[j] = i1
                
                dp1 = norms[i]
                norms[i] = norms[j]
                norms [j] = dp1

    for i in xrange (0, cmsize):
        for j in xrange (i, cmsize) :
            cm [i][j] = cm2[indexes[i]][indexes[j]]
            cm[j][i] = cm[i][j]
                                                                                                                                                                                                                  
    return cm
