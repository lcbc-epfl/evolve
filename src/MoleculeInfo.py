'''
main.py

@author: Nicholas Browning
'''

import openbabel
from openbabel import OBResidue
from openbabel import OBAtom
import numpy as np
import constants
from Cython.Compiler.Nodes import ContinueStatNode

# defines the expect sidechain atoms for the aminoacids (copy&paste form Basilisk CalcAngles)
sidechain = {
         "ALA": ['cb'],
         "ARG": ['cb', 'cg', 'cd', 'ne', 'cz', 'nh1', 'nh2'],
         "ASN": ['cb', 'cg', 'od1', 'nd2'],
         "ASP": ['cb', 'cg', 'od1', 'od2'],
         "CYS": ['cb', 'sg'],
         "GLU": ['cb', 'cg', 'cd', 'oe1', 'oe2'],
         "GLN": ['cb', 'cg', 'cd', 'oe1', 'ne2'],
         "GLY": [],
         "HIS": ['cb', 'cg', 'nd1', 'cd2', 'ce1', 'ne2'],
         "ILE": ['cb', 'cg1', 'cg2', 'cd1'],
         "LEU": ['cb', 'cg', 'cd1', 'cd2'],
         "LYS": ['cb', 'cg', 'cd', 'ce', 'nz'],
         "MET": ['cb', 'cg', 'sd', 'ce'],
         "PHE": ['cb', 'cg', 'cd1', 'cd2', 'ce1', 'ce2', 'cz'],
         "PRO": ['cb', 'cg', 'cd'],
         "SER": ['cb', 'og'],
         "THR": ['cb', 'og1', 'cg2'],
         "TRP": ['cb', 'cg', 'cd1', 'cd2', 'ne1', 'ce2', 'ce3', 'cz2', 'cz3', 'ch2'],
         "TYR": ['cb', 'cg', 'cd1', 'cd2', 'ce1', 'ce2', 'cz', 'oh'],
         "VAL": ['cb', 'cg1', 'cg2']
         }

# defines the different atoms needed to calculate the chi angles
# for all the different residue types ... well at least the common 
# 20 ones (copy&paste form Basilisk CalcAngles)
chi_atoms = {  
         "ALA": {},
         "ARG": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd"],
                 "x3": ["cb", "cg", "cd", "ne"], "x4": ["cg", "cd", "ne", "cz"]},
         "ASN": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "od1"]},
         "ASP": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "od1"]},
         "CYS": {"x1": ["n", "ca", "cb", "sg"]},
         "GLU": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd"],
                 "x3": ["cb", "cg", "cd", "oe1"]},
         "GLN": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd"],
                 "x3": ["cb", "cg", "cd", "oe1"]},
         "GLY": {},
         "HIS": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd2"]},
         "ILE": {"x1": ["n", "ca", "cb", "cg1"], "x2": ["ca", "cb", "cg1", "cd1"]},
         "LEU": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd1"]},
         "LYS": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd"],
                 "x3": ["cb", "cg", "cd", "ce"], "x4": ["cg", "cd", "ce", "nz"]},
         "MET": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "sd"],
                 "x3": ["cb", "cg", "sd", "ce"]},
         "PHE": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd1"]},
         "PRO": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd"],
                 "x3": ["cb", "cg", "cd", "n"], "x4": ["cg", "cd", "n", "ca"]},
         "SER": {"x1": ["n", "ca", "cb", "og"]},
         "THR": {"x1": ["n", "ca", "cb", "og1"]},
         "TRP": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd1"]},
         "TYR": {"x1": ["n", "ca", "cb", "cg"], "x2": ["ca", "cb", "cg", "cd1"]},
         "VAL": {"x1": ["n", "ca", "cb", "cg1"]}
         }

def get_chi_atoms(obres):
    """
    For each residuetype there is a certain set of dihedral angles in the sidechain 
    which again are defined by a certain set of atoms each.
    
    @return: returns a dictionary of arrays, describing the atoms for all the angles 
    @rtype: dictionary
    """
    residuetype = obres.GetName()
    if residuetype.upper() not in chi_atoms:
        sys.stderr.write("Warning: Unknown residuetype " + residuetype + " in residue list in MoleculeInfo.py\n")
        return {}
    #
    out = chi_atoms[residuetype.upper()]
    for k in range(0, len(out)):
        out[k] = out[k].upper()
    
    return out

def get_sidechain_atoms_names(obres):
    residuetype = obres.GetName()
    if residuetype.upper() not in sidechain:
        sys.stderr.write("Warning: Unknown residuetype " + residuetype + " in residue list in MoleculeInfo.py\n")
        return {}
    #
    out = sidechain[residuetype.upper()]
    for k in range(0, len(out)):
        out[k] = out[k].upper()
    
    return out

def get_sidechain_atoms(obres):
    # return a list of three lists, 
    # first list is the atoms' names, second is the corresponding OBatoms, 
    # third is the corresponding (idx) indices in the overall molecule (OBmol)
    OBatoms_in_res = []
    IDatoms_in_res = []
    for obatom in openbabel.OBResidueAtomIter(obres):
        OBatoms_in_res.append(obatom)
        IDatoms_in_res.append(obres.GetAtomID(obatom).strip())

    sidechain_IDatoms = get_sidechain_atoms_names(obres)
    sidechain_OBatoms = []
    sidechain_NUMatoms = []
    for k in range(0, len(sidechain_IDatoms)):
        atomID = sidechain_IDatoms[k]
        loc = IDatoms_in_res.index(atomID)
        sidechain_OBatoms.append(OBatoms_in_res[loc])
        sidechain_NUMatoms.append(OBatoms_in_res[loc].GetIdx())

    return (sidechain_IDatoms, sidechain_OBatoms, sidechain_NUMatoms)

def getResType(obres):
    '''get 3 letter code of OBRes - needed to parse richardson library'''
    res = obres.GetName()
    
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
    
    # need to do additional tests for Glycine
    res = getResType(obres)
    for obatom in openbabel.OBAtomAtomIter(alpha_carbon):
        if (res == "GLY"):
            if (obatom.GetType() == 'H'):
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
        # if alpha_carbon != None:
            # print('CA: {}'.format(alpha_carbon.GetId()))
        
        bb_nitrogen = getBBNitrogen(obres)
        # if bb_nitrogen != None:
            # print('N: {}'.format(bb_nitrogen.GetId()))
        
        neg1_carboxl = getNeg1BBCarboxyl(obres)
        # if neg1_carboxl != None:
            # print('C-1: {}'.format(neg1_carboxl.GetId()))

        carboxl = getBBCarboxyl(obres)
        # if carboxl != None:
            # print('C: {}'.format(carboxl.GetId()))
        
        plus1_nitrogen = getPlus1BBNitrogen(obres)
        # if plus1_nitrogen != None:
            # print('N+1: {}'.format(plus1_nitrogen.GetId()))
        
        
        phi = 999
        psi = 999
        
        if (neg1_carboxl is not None and carboxl is not None):
            phi = mol.GetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl)
        if (bb_nitrogen is not  None and plus1_nitrogen is not None):
            psi = mol.GetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen)
        
        result.append((phi, psi))
    return result     


def setChi1DihedralAngle(mol, residue_index, angle_deg):
    
    obres = mol.GetResidue(residue_index)
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    
    gamma_atom = getChi1DihedralAtom(obres)
    
    if (gamma_atom is None):
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
    
    if (gamma_atom is None):
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
    
    if (neg1_carboxl is not None and carboxl and phi != 999):

        mol.SetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl, phi * (np.pi / 180.0))
    
    if (bb_nitrogen is not None and plus1_nitrogen and psi != 999):
    
        mol.SetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen, psi * (np.pi / 180.0))

     
def getPlus1BBNitrogen(obres):
    
    bb_n = getBBCarboxyl(obres)
    
    for obatom in openbabel.OBAtomAtomIter(bb_n):

        if (obatom.IsNitrogen()):
            return obatom


def get_atoms_per_residue(mol):
    '''
    Gives a list of lists of j atom character belonging to residue i as well as the OBAtoms equivalent
    '''
    sorted_atoms_ID = [[] for _ in xrange(mol.NumResidues())]
    sorted_atoms_OB = [[] for _ in xrange(mol.NumResidues())]

    for i, ob_res in enumerate(openbabel.OBResidueIter(mol)):
        # print "Res", i, ob_res.GetName()
        
        for ob_atom in openbabel.OBResidueAtomIter(ob_res):
            sorted_atoms_OB[i].append(ob_atom)
            sorted_atoms_ID[i].append(ob_res.GetAtomID(ob_atom).strip())

        # print(sorted_atoms_ID[i])

    return (sorted_atoms_ID, sorted_atoms_OB)



if __name__ == '__main__':
    import sys
    print "Starting"

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, sys.argv[1]) 
    
    
    for j in xrange(1, mol.NumAtoms() + 1): 
        print debugAtom(mol.GetAtom(j))
        
        
    for j in xrange(0, mol.NumResidues()): 
        
        obres = mol.GetResidue(j)
        
        print "--", obres.GetName(), "--"
        
        print get_sidechain_atoms(obres)