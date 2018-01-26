'''
main.py

@author: Nicholas Browning
'''

import openbabel
from openbabel import OBResidue
from openbabel import OBAtom
import numpy as np
import constants


def getResType(obres):
    '''get 3 letter code of OBRes - needed to parse richardson library'''
    res = obres.GetName()
    
    res = res.strip()
    res = res[:3]
    # check richardson rotamers and return rotamer_type
    for i in xrange (0, len(constants.rotamers)):
        if (constants.rotamers[i] == res):
            return constants.rotamer_types[i]
        
    # seems to be a standard AA, search 3-letter code list
    for i in xrange (0, len(constants.rotamer_types)):
        if (constants.rotamer_types[i] == res):
            return constants.rotamer_types[i]
    
    if (res == "PRO"):
        return "PRO"
    # non standard residue?
    return None

def getAlphaCarbon(obres):
    for obatom in openbabel.OBResidueAtomIter(obres):
        # atomID = obres.GetAtomID(obatom)
        
        # if (atomID.strip(" ") == "CA"):
        #    return obatom
        
        if (obatom.IsCarbon()):
           
            
            carboxl_carbon = getConnectedCarboxylCarbon(obatom)
            # can identify carboxyl carbon in any peptide
            if (carboxl_carbon is not None):
                # now can find proline nitrogen
                amide_nitro = getConnectedAmideNitrogen(obatom)
                
                if (amide_nitro is not None):
                    return obatom
    return None
        
def getConnectedAmideNitrogen(carbon_atom):
    
    obres = carbon_atom.GetResidue()
    # can have trouble here for proline due to it's connectivity, need separate rules
    #    --> assume that atom is the alpha_carbon
    if (obres.GetName() == "PRO"):
        for obatom in openbabel.OBAtomAtomIter(carbon_atom):
            if (obatom.IsNitrogen()):

                return obatom
                    
    else:
        for obatom in openbabel.OBAtomAtomIter(carbon_atom):
            if (obatom.IsNitrogen()):
              
                for obatom2 in openbabel.OBAtomAtomIter(obatom):
                    if (obatom2.IsHydrogen()):
                     
                        return obatom
    return None

def getConnectedCarboxylCarbon(atom):
    
    obres = atom.GetResidue()
    
    for obatom in openbabel.OBAtomAtomIter(atom):
        if (obatom.IsCarbon()):
       
            for obbond in openbabel.OBAtomBondIter(obatom):
                if (countBonds(obbond.GetNbrAtom(obatom)) == 1 and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                    return obatom
    return None
        
        
def getBetaAtom(obres):
    alpha_carbon = getAlphaCarbon(obres)

    bbcarboxl = getBBCarboxyl(obres)
 
    if alpha_carbon is None:
        return None
    
    res = getResType(obres)
    for obatom in openbabel.OBAtomAtomIter(alpha_carbon):
        
        if (res == "GLY"):

            bbnitrogen = getBBNitrogen(obres)
    
            if (obatom.GetAtomicNum() == 1):
                
                # select the L enantiomer hydrogen
                
                h_vec = np.asarray([obatom.GetX() - bbcarboxl.GetX(), obatom.GetY() - bbcarboxl.GetY(), obatom.GetZ() - bbcarboxl.GetZ()])
                ca_vec = np.asarray([alpha_carbon.GetX() - bbcarboxl.GetX(), alpha_carbon.GetY() - bbcarboxl.GetY(), alpha_carbon.GetZ() - bbcarboxl.GetZ()])
                n_vec = np.asarray([bbnitrogen.GetX() - alpha_carbon.GetX(), bbnitrogen.GetY() - alpha_carbon.GetY(), bbnitrogen.GetZ() - alpha_carbon.GetZ()])

                orth_vec = np.cross(ca_vec, n_vec)
        
                if (angle(h_vec, orth_vec) < 90.0):
                    # this corresponds to the L enantiomer
                    # print angle(h_vec, orth_vec)
                    return obatom
                

        else:
            if (obatom.IsCarbon() and obatom.GetIdx() != bbcarboxl.GetIdx()):
                return obatom
    return None
        
def getBBNitrogen(obres):
    alpha_carbon = getAlphaCarbon(obres)
    
    if alpha_carbon is None:
        return None
    
    for obatom in openbabel.OBAtomAtomIter(alpha_carbon):
        if (obatom.IsNitrogen()):
            return obatom
    return None
        
def getNeg1BBCarboxyl(obres):
    ca = getAlphaCarbon(obres)
    bb_n = getBBNitrogen(obres)
    
    for obatom in openbabel.OBAtomAtomIter(bb_n):
        if (obatom == ca):
            continue
        
        if (obatom.IsHydrogen()):
            continue
        
        for obbond in openbabel.OBAtomBondIter(obatom):
            if (countBonds(obbond.GetNbrAtom(obatom)) == 1  and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                return obatom
    return None
            
def getBBCarboxyl(obres):
    ca_atom = getAlphaCarbon(obres)
    
    # print debugAtom(ca_atom)
    for obatom in openbabel.OBAtomAtomIter(ca_atom):
        
        # print obatom.GetType()

        if (not obatom.IsCarbon()):
            continue
        
        # print "--", obatom.GetType()
        
        for obbond in openbabel.OBAtomBondIter(obatom):
            # print obbond.GetBO(), obbond.GetBeginAtom().GetType(), obbond.GetEndAtom().GetType()
            if (countBonds(obbond.GetNbrAtom(obatom)) == 1  and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                return obatom
    return None
            

def getPhiPsiDihedrals(mol, residue_indexes):
    result = []
    
    for i in xrange (0, len(residue_indexes)):
        
        obres = mol.GetResidue(residue_indexes[i])
        
        alpha_carbon = getAlphaCarbon(obres)
        bb_nitrogen = getBBNitrogen(obres)
        neg1_carboxl = getNeg1BBCarboxyl(obres)
        carboxl = getBBCarboxyl(obres)
        plus1_nitrogen = getPlus1BBNitrogen(obres)

        phi = 999
        psi = 999
        
        if (neg1_carboxl != None and carboxl != None):
            phi = mol.GetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl)
        if (bb_nitrogen != None and plus1_nitrogen != None):
            psi = mol.GetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen)
        
        result.append((phi, psi))
    return result     


def setChi1DihedralAngle(mol, residue_index, angle_deg):
    
    obres = mol.GetResidue(residue_index)
    
    obres = mol.GetResidue(residue_index)
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    
    gamma_atom = getChi1DihedralAtom(obres)
    
    if (gamma_atom == None):
        return False
    
    mol.SetTorsion(bb_nitrogen, alpha_carbon, beta_atom, gamma_atom, angle_deg * (np.pi / 180.0))

    return True

def getPosition(obatom):
    return np.asarray([obatom.GetX(), obatom.GetY(), obatom.GetZ()])

def getChi1DihedralAngle(mol, obres):
    
    
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    gamma_atom = getChi1DihedralAtom(obres)
 
    # print debugAtom(alpha_carbon)
    # print debugAtom(bb_nitrogen)
    # print debugAtom(beta_atom)
    # print debugAtom(gamma_atom)
    
    if (gamma_atom == None):
        return None
    
    return mol.GetTorsion(bb_nitrogen, alpha_carbon, beta_atom, gamma_atom)

def countBonds(obatom):
    count = 0
    for i in xrange(1, 4):
        count += obatom.CountBondsOfOrder(i)
    return count

def debugAtom(obatom):
    return "AN:", obatom.GetAtomicNum(), "T:", obatom.GetType(), "idx:", obatom.GetIdx(), "id:", obatom.GetId(), "BO(1):", obatom.CountBondsOfOrder(1), obatom.GetVector()

def debugBond(obbond):
    return "begin_idx:", obbond.GetBeginAtomIdx(), "end_idx:", obbond.GetEndAtomIdx(), "norm:", obbond.GetLength(), "BO:", obbond.GetBO()

  

def getChi1DihedralAtom(obres):
    
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    gamma_atom = None
    
    if beta_atom is None:
        return None
    
    res = getResType(obres)
    
    
    for obatom in openbabel.OBAtomAtomIter(beta_atom):
    
        if (obatom.GetIdx() == alpha_carbon.GetIdx()):
            continue
        
        if ((res == "SER" or res == "THR") and obatom.IsOxygen()):
            return obatom
        elif (res == "CYS" and obatom.IsSulphur()):
            return obatom
        elif (res == "ALA" and obatom.IsHydrogen()):
            return obatom
        elif((res != "SER" and res != "THR" and res != "CYS" and res != "ALA") and obatom.IsCarbon()):
            return obatom
        
    return None
    
            
def getAllPhiPsiDihedrals(mol):
    return getPhiPsiDihedrals(mol, np.arange(mol.NumResidues()))

def SetPhiPsi(mol, obres, phi, psi):
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    neg1_carboxl = getNeg1BBCarboxyl(obres)
    carboxl = getBBCarboxyl(obres)
    plus1_nitrogen = getPlus1BBNitrogen(obres)
    
    if (neg1_carboxl and carboxl and phi != 999):
        mol.SetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl, phi * (np.pi / 180.0))
    
    if (bb_nitrogen and plus1_nitrogen and psi != 999):
        mol.SetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen, psi * (np.pi / 180.0))
     
def getPlus1BBNitrogen(obres):
    
    bb_n = getBBCarboxyl(obres)
    
    for obatom in openbabel.OBAtomAtomIter(bb_n):

        if (obatom.IsNitrogen()):
            return obatom



# ATOM    148  N   GLY     9      38.893  23.411  33.937  0.00  0.00
# ATOM    150  CA  GLY     9      39.937  23.667  34.946  0.00  0.00
# ATOM    151  HA2 GLY     9      40.563  24.496  34.619  0.00  0.00
# ATOM    152  HA3 GLY     9      39.525  24.058  35.876  0.00  0.00
# ATOM    153  C   GLY     9      40.820  22.504  35.210  0.00  0.00

n_vec = np.asarray([38.893  , 23.411, 33.937])
c_vec = np.asarray([40.820 , 22.504 , 35.210])

ca_vec = np.asarray([39.937, 23.667, 34.946])
ha2_vec = np.asarray([40.563, 24.496, 34.619])
ha3_vec = np.asarray([39.525 , 24.058 , 35.876])


def angle(vec_1, vec_2):
    return np.arccos(np.dot(vec_1, vec_2) / (np.linalg.norm(vec_1) * np.linalg.norm(vec_2))) * 180 / np.pi


if __name__ == '__main__':
#     print "Starting"
# 
#     obConversion = openbabel.OBConversion()
#     obConversion.SetInAndOutFormats("pdb", "pdb")
# 
#     mol = openbabel.OBMol()
#     obConversion.ReadFile(mol, "gpgg.pdb") 
#     
#     print mol.NumAtoms()
#     print mol.NumBonds()
#     print mol.NumResidues()
#     
#     # getAllPhiPsiDihedrals(mol)
#         
#     for j in xrange(0, mol.NumResidues()): 
#         
#         obres = mol.GetResidue(j)
#         
#         print "-----"
#         for obatom in openbabel.OBResidueAtomIter(obres):
#             print obatom.GetIdx(), obatom.GetAtomicNum()
#         print "-----"   
#         alpha_carbon = getAlphaCarbon(obres)
#         
#         print "ALpha:", alpha_carbon
#         
#         bb_nitrogen = getBBNitrogen(obres)
#         neg1_carboxl = getNeg1BBCarboxyl(obres)
#         carboxl = getBBCarboxyl(obres)
#         plus1_nitrogen = getPlus1BBNitrogen(obres)
#         
#         if (not neg1_carboxl):
#             continue
#         
#         print j, obres.GetName(), "::" , neg1_carboxl, "---", bb_nitrogen, "---", alpha_carbon, "---", carboxl, "---", plus1_nitrogen
#         
#         mol.SetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl, -60 * (np.pi / 180.0))
#         
#     print "Writing to file"
#     angles = getAllPhiPsiDihedrals(mol)
#     print angles
#     obConversion.WriteFile(mol, "test.pdb")

    
    h2_vec = ha2_vec - ca_vec
    h3_vec = ha3_vec - ca_vec
    
    orth_vec = np.cross(ca_vec - c_vec, n_vec - ca_vec)
    
    
    print angle(h2_vec, c_vec - ca_vec)
    print angle(h3_vec, c_vec - ca_vec)
    
    print angle(h2_vec, orth_vec)
    print angle(h3_vec, orth_vec)
    
