'''
MoleculeInfo.py

Contains helper functions for finding + identifying certain atoms in proteins.

.. codeauthor:: Nicholas Browning
'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import str
from builtins import range
from past.utils import old_div
# import openbabel
# from openbabel import OBResidue
# from openbabel import OBAtom
from openbabel import openbabel
from openbabel.openbabel import OBResidue
from openbabel.openbabel import OBAtom
import numpy as np
from . import constants

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
# 20 ones (copy&paste from Basilisk CalcAngles)
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
         "VAL": {"x1": ["n", "ca", "cb", "cg1"]}}


def getResType(obres):
    """
    get 3 letter code of an OBResidue object
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    res = obres.GetName()

    # check richardson rotamers and return rotamer_type
    for i in range (0, len(constants.rotamers)):
        if (constants.rotamers[i] == res):
            return constants.rotamer_types[i]
        
    # seems to be a standard AA, search 3-letter code list
    # DUPLICATE
    for i in range (0, len(constants.rotamer_types)):
        if (constants.rotamer_types[i] == res):
            return constants.rotamer_types[i]
    
    if (res == "PRO"):
        return "PRO"
    # non standard residue?
    return None


def getAlphaCarbon(obres):
    """
    Finds the alpha carbon in a given residue.
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    res = obres.GetName()
    for obatom in openbabel.OBResidueAtomIter(obres):

        if (obatom.GetAtomicNum()==6):
            carboxl_carbon = getConnectedCarboxylCarbon(obatom)
            # identify carboxyl carbon in peptide
            if (carboxl_carbon is not None):
                # now also find the nitrogen
                amide_nitro = getConnectedAmideNitrogen(obatom)
                if (amide_nitro is not None):
                    return obatom
    return None

        
def getConnectedAmideNitrogen(calpha_atom):
    """
    returns the OBAtom corresponding to the amide nitrogen along a protein backbone. Similar to getBBNitrogen but with additional checks.
    
    Parameters
    ----------
    calpha_atom : OBAtom
        OBAtom corresponding to the alpha carbon
    """
    
    obres = calpha_atom.GetResidue()
    # can have trouble here for proline due to it's connectivity, need separate rules
    #    --> assume that atom is the alpha_carbon
    if (obres.GetName() == "PRO"):
        for obatom in openbabel.OBAtomAtomIter(calpha_atom):
            if (obatom.GetAtomicNum()==7):
                return obatom
                    
    else:
        for obatom in openbabel.OBAtomAtomIter(calpha_atom):
            if (obatom.GetAtomicNum()==7):
                for obatom2 in openbabel.OBAtomAtomIter(obatom):
                    if (obatom2.GetAtomicNum()==1):
                        return obatom
    return None


def getConnectedCarboxylCarbon(atom):
    """
    returns the OBAtom corresponding to the carboxyl carbon along a protein backbone. Similar to getBBCarboxyl but with additional checks.
    
    Parameters
    ----------
    calpha_atom : OBAtom
        OBAtom corresponding to the alpha carbon
    """
    obres = atom.GetResidue()
    
    for obatom in openbabel.OBAtomAtomIter(atom):
        if (obatom.GetAtomicNum()==6):
       
            for obbond in openbabel.OBAtomBondIter(obatom):
                if (countBonds(obbond.GetNbrAtom(obatom)) == 1 and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                    return obatom
    return None
        
        
def getBetaAtom(obres):
    """
     Returns the atom connected to the alpha carbon. For GLY, this corresponds to the hydrogen leading to the L-enantiomer, otherwise for non-GLY residues
     it corresponds to the CB atom.
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """

    alpha_carbon = getAlphaCarbon(obres)
    
    bbcarboxyl = getBBCarboxyl(obres)
    
    bbnitrogen = getBBNitrogen(obres)
    
    if alpha_carbon is None:
        return None
    
    if bbcarboxyl is None:
        return None
    
    if bbnitrogen is None:
        return None
    
    alpha_vec = np.asarray([alpha_carbon.GetX(), alpha_carbon.GetY(), alpha_carbon.GetZ()])
    nitro_vec = np.asarray([bbnitrogen.GetX(), bbnitrogen.GetY(), bbnitrogen.GetZ()])
    carboxyl_vec = np.asarray([bbcarboxyl.GetX(), bbcarboxyl.GetY(), bbcarboxyl.GetZ()])
    
    res = getResType(obres)
    for obatom in openbabel.OBAtomAtomIter(alpha_carbon):
        if (res == "GLY"):
            if (obatom.GetType() == 'H'):
                 # want to find the L-enantiomer of the two GLY hydrogens on CA
                ob_vec = np.asarray([obatom.GetX(), obatom.GetY(), obatom.GetZ()])
                
                RH = ob_vec - alpha_vec
                RC = carboxyl_vec - alpha_vec
                RN = nitro_vec - alpha_vec
                perp_vec = np.cross(RN, RC)
               
                if (np.dot(RH, perp_vec) > 0.0):
                    return obatom

        else:
            if (bbcarboxyl is not None and obatom.GetAtomicNum()==6 and obatom.GetIdx() != bbcarboxyl.GetIdx()):
                return obatom
    return None

        
def getBBNitrogen(obres):
    """
     Returns OBAtom corresponding to the backbone nitrogen in a given OBResidue
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    alpha_carbon = getAlphaCarbon(obres)
    
    if (alpha_carbon is None):
        return None
    
    for obatom in openbabel.OBAtomAtomIter(alpha_carbon):
        if (obatom.GetAtomicNum()==7):
            return obatom
    return None

        
def getNeg1BBCarboxyl(obres):
    """
     Returns the OBAtom corresponding to the backbone nitrogen of the preceding residue connected to the given OBResidue

    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    ca = getAlphaCarbon(obres)
    bb_n = getBBNitrogen(obres)
    
    if (ca is None and bb_n is None):
        return None
    
    for obatom in openbabel.OBAtomAtomIter(bb_n):
        if (obatom == ca):
            continue
        
        if (obatom.GetAtomicNum()==1):
            continue
        
        for obbond in openbabel.OBAtomBondIter(obatom):
            if (countBonds(obbond.GetNbrAtom(obatom)) == 1  and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                return obatom
    return None

            
def getBBCarboxyl(obres):
    """
     Returns the OBAtom corresponding to the backbone carboxyl carbon of a given OBResidue
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    ca_atom = getAlphaCarbon(obres)
    
    if (ca_atom is None):
        return None
    
    for obatom in openbabel.OBAtomAtomIter(ca_atom):
        
        if (not obatom.GetAtomicNum()==6):
            continue
        
        for obbond in openbabel.OBAtomBondIter(obatom):
            if (countBonds(obbond.GetNbrAtom(obatom)) == 1  and obbond.GetNbrAtom(obatom).GetAtomicNum() == 8):
                return obatom
    return None


def getPlus1BBNitrogen(obres):
    """
    Returns the OBAtom corresponding to the backbone nitrogen of the next residue connected to the given OBResidue
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    bb_n = getBBCarboxyl(obres)
    
    if (bb_n is None):
        return None
    
    for obatom in openbabel.OBAtomAtomIter(bb_n):

        if (obatom.GetAtomicNum()==7):
            return obatom
    return None


def getChiDihedralsByAtomIndex(mol, chi_dihedral_atom_indexes):
    
    """
    Returns the set of chi-1 angles for a given input of [N, 4] atom indexes, where N is the number of query residues
    
    Parameters
    ----------
    mol : OBMol
        protein
    chi_dihedral_atom_indexes: [N, 4] matrix
        matrix with N query residues, where each column corresponds to the atom indexes belonging to the dihedral
    """
        
    chis = []
    if (chi_dihedral_atom_indexes is None or len(chi_dihedral_atom_indexes) == 0):
        return chis

    for j in range(len(chi_dihedral_atom_indexes)):
        
        res_chis = []
        
        if len(chi_dihedral_atom_indexes[j]) == 0:
            chis.append(res_chis)
            continue

        for k in range(0, len(chi_dihedral_atom_indexes[j])):
        
            chi_angle_atoms = chi_dihedral_atom_indexes[j][k]
            res_chis.append(mol.GetTorsion(chi_angle_atoms[0], chi_angle_atoms[1], chi_angle_atoms[2], chi_angle_atoms[3]))
            
        chis.append(res_chis)
        
    return np.asarray(chis)

                
def getChiDihedralAtomIndexes(mol, residue_indexes):
    
    """
    Returns a list of all atoms involved in sidechain chi-1 dihedrals for a set of residues
    
    Parameters
    ----------
    mol : OBMol
        protein
    residue_indexes: list (int)
        list of residue indexes
    """

    chi_atoms = []
    
    for i in range (0, len(residue_indexes)):
        
        obres = mol.GetResidue(residue_indexes[i])
            
        obres = mol.GetResidue(i)
        sidechain_atoms_dict = get_chi_atoms(obres)
        
        if (not sidechain_atoms_dict or obres.GetName() == "PRO"):
            chi_atoms.append([])
            continue
            
        num_chi_angles = len(sidechain_atoms_dict)

        res_chi_atoms = []
            
        for j in range(num_chi_angles):
            atoms = [None, None, None, None]
        
            sidechain_torsion_atomsID = [atomID.upper() for atomID in sidechain_atoms_dict['x' + str(j + 1)]]
                
            (IDatoms_in_res, OBatoms_in_res, NUMatoms_in_res) = get_atoms(obres)
                
            for k in range(0, 4):
                atoms[k] = OBatoms_in_res[IDatoms_in_res.index(sidechain_torsion_atomsID[k])].GetIdx()
        
            res_chi_atoms.append(atoms)
                
        chi_atoms.append(res_chi_atoms)
            
    return np.asarray(chi_atoms)


def getPhiPsiDihedralByAtomIndex(mol, dihedral_atom_indexes):
    """
    Returns a list of [phi, psi] dihedrals for the backbone dihedrals given by dihedral_atom_indexes
    
    Parameters
    ----------
    mol : OBMol
        protein
    dihedral_atom_indexes: [N,2,4] list 
        list of atom indexes involved in both the phi and psi dihedrals
    """
    
    result = []
    
    for i in range(len(dihedral_atom_indexes)):
       
        phi = None
        psi = None
        
        if len(dihedral_atom_indexes[i][0]) != 0:
            phi = mol.GetTorsion(dihedral_atom_indexes[i][0][0], dihedral_atom_indexes[i][0][1], dihedral_atom_indexes[i][0][2], dihedral_atom_indexes[i][0][3])
        
        if len(dihedral_atom_indexes[i][1]) != 0:
            psi = mol.GetTorsion(dihedral_atom_indexes[i][1][0], dihedral_atom_indexes[i][1][1], dihedral_atom_indexes[i][1][2], dihedral_atom_indexes[i][1][3])  
        
        result.append([phi, psi])
        
    return np.asarray(result) 


def getPhiPsiDihedralAtomIndexes(mol, residue_indexes): 
    """
    Returns a list of all atoms involved in backbone phi, psi dihedrals for a set of residues
    
    Parameters
    ----------
    mol : OBMol
        protein
    residue_indexes: list (int)
        list of residue indexes
    """
    
    result = []
    
    for i in range (0, len(residue_indexes)):
        
        dihedrals = []
        
        obres = mol.GetResidue(residue_indexes[i])
        
        alpha_carbon = getAlphaCarbon(obres)
        
        bb_nitrogen = getBBNitrogen(obres)
      
        neg1_carboxl = getNeg1BBCarboxyl(obres)
    
        carboxl = getBBCarboxyl(obres)
     
        plus1_nitrogen = getPlus1BBNitrogen(obres)
       
        if (neg1_carboxl is not None and bb_nitrogen is not None and alpha_carbon is not None and carboxl is not None):
            dihedrals.append([neg1_carboxl.GetIdx(), bb_nitrogen.GetIdx(), alpha_carbon.GetIdx(), carboxl.GetIdx()])
        else:
            dihedrals.append([])
            
        if (bb_nitrogen is not None and alpha_carbon is not None and carboxl is not None and plus1_nitrogen is not None):
            dihedrals.append([bb_nitrogen.GetIdx(), alpha_carbon.GetIdx(), carboxl.GetIdx(), plus1_nitrogen.GetIdx()])
        else:
            dihedrals.append([])
        
        result.append(dihedrals)

    return np.asarray(result)


def getPhiPsiDihedrals(mol, residue_indexes):
    """
    Returns a [N,2] list of all backbone phi, psi dihedrals for a set of residues
    
    Parameters
    ----------
    mol : OBMol
        protein
    residue_indexes: list (int)
        list of residue indexes
    """
    
    result = []

    for i in range (0, len(residue_indexes)):
      
        obres = mol.GetResidue(residue_indexes[i])

        alpha_carbon = getAlphaCarbon(obres)
        bb_nitrogen = getBBNitrogen(obres)
        neg1_carboxl = getNeg1BBCarboxyl(obres)
        carboxl = getBBCarboxyl(obres)
        plus1_nitrogen = getPlus1BBNitrogen(obres)
        
        phi = 999
        psi = 999
        
        if (neg1_carboxl is not None and carboxl is not None):
            phi = mol.GetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl)
        if (bb_nitrogen is not  None and plus1_nitrogen is not None):
            psi = mol.GetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen)
        
        result.append((phi, psi))
    return result     


def setChi1DihedralAngle(mol, residue_index, angle_deg):
    """
    Sets the first chi dihedral (in degrees) for a given residue
    
    Parameters
    ----------
    mol : OBMol
        protein
    residue_index: int
        residue index
    angle_deg:
        angle in degrees
    """
    
    obres = mol.GetResidue(residue_index)
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    
    gamma_atom = getChi1DihedralAtom(obres)
    
    if (gamma_atom is None):
        return False
    
    mol.SetTorsion(bb_nitrogen, alpha_carbon, beta_atom, gamma_atom, angle_deg * (old_div(np.pi, 180.0)))

    return True


def getPosition(obatom):
    return np.asarray([obatom.GetX(), obatom.GetY(), obatom.GetZ()])


def getChi1DihedralAngle(mol, obres):
    """
    Gets the first chi dihedral (in degrees) for a given OBResidue
    
    Parameters
    ----------
    mol : OBMol
        protein
    obres: OBResidue
        query residue
    """
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    beta_atom = getBetaAtom(obres)
    gamma_atom = getChi1DihedralAtom(obres)
    
    if (gamma_atom is None):
        return None
    
    return mol.GetTorsion(bb_nitrogen, alpha_carbon, beta_atom, gamma_atom)


def countBonds(obatom):
    count = 0
    for i in range(1, 4):
        count += obatom.CountBondsOfOrder(i)
    return count


def debugAtom(obatom):
    return "AN:", obatom.GetAtomicNum(), "T:", obatom.GetType(), "idx:", obatom.GetIdx(), "id:", obatom.GetId(), "BO(1):", obatom.CountBondsOfOrder(1), obatom.GetVector()


def debugBond(obbond):
    return "begin_idx:", obbond.GetBeginAtomIdx(), "end_idx:", obbond.GetEndAtomIdx(), "norm:", obbond.GetLength(), "BO:", obbond.GetBO()
  

def getChi1DihedralAtom(obres):
    """
    Gets the final atom involved in the chi-1 dihedral for a given OBResidue
    
    Parameters
    ----------
    mol : OBMol
        protein
    obres: OBResidue
        query residue
    """
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
        
        if ((res == "SER" or res == "THR") and obatom.GetAtomicNum()==8):
            return obatom
        elif (res == "CYS" and obatom.GetAtomicNum()==16):
            return obatom
        elif (res == "ALA" and obatom.GetAtomicNum()==1):
            return obatom
        elif((res != "SER" and res != "THR" and res != "CYS" and res != "ALA") and obatom.GetAtomicNum()==6):
            return obatom
        
    return None
    
            
def getAllPhiPsiDihedrals(mol):
    """
    Returns a list of all [phi, psi] dihedrals in a protein
    
    Parameters
    ----------
    mol : OBMol
        protein
    """
    return getPhiPsiDihedrals(mol, np.arange(mol.NumResidues()))


def SetPhiPsi(mol, obres, phi, psi):
    """
    Sets a specific [phi, psi] angle (in degrees) for a given residue
    
    Parameters
    ----------
    mol : OBMol
        protein
    obres: OBResidue
        query residue
    phi: float
        phi angle in degrees
    psi: float
        psi angle in degrees
    """
    
    alpha_carbon = getAlphaCarbon(obres)
    bb_nitrogen = getBBNitrogen(obres)
    neg1_carboxl = getNeg1BBCarboxyl(obres)
    carboxl = getBBCarboxyl(obres)
    plus1_nitrogen = getPlus1BBNitrogen(obres)
    
    if (neg1_carboxl is not None and carboxl and phi != 999):

        mol.SetTorsion(neg1_carboxl, bb_nitrogen, alpha_carbon, carboxl, phi * (old_div(np.pi, 180.0)))
    
    if (bb_nitrogen is not None and plus1_nitrogen and psi != 999):
    
        mol.SetTorsion(bb_nitrogen, alpha_carbon, carboxl, plus1_nitrogen, psi * (old_div(np.pi, 180.0)))


def get_atoms_per_residue(mol):
    """
    Returns a tuple of all atom IDs and their OBAtoms objects beloging to each residue in a given protein
    
    Parameters
    ----------
    mol : OBMol
        protein
    """
    
    sorted_atoms_ID = [[] for _ in range(mol.NumResidues())]
    sorted_atoms_OB = [[] for _ in range(mol.NumResidues())]

    for i, ob_res in enumerate(openbabel.OBResidueIter(mol)):

        for ob_atom in openbabel.OBResidueAtomIter(ob_res):
            sorted_atoms_OB[i].append(ob_atom)
            sorted_atoms_ID[i].append(ob_res.GetAtomID(ob_atom).strip())

    return (sorted_atoms_ID, sorted_atoms_OB)

# defines the expect sidechain atoms for the aminoacids (copy&paste form Basilisk CalcAngles)


def get_sidechain_atoms_names(obres):
    """
    TODO
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    residuetype = obres.GetName()
    if residuetype.upper() not in sidechain:
        sys.stderr.write("Warning: Unknown residuetype " + residuetype + " in residue list in MoleculeInfo.py\n")
        return {}
    #
    out = sidechain[residuetype.upper()]
    for k in range(0, len(out)):
        out[k] = out[k].upper()
    
    return out


def get_atoms(obres):
    """
    TODO
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    # return a list of three lists, 
    # first list is the atoms' names, second is the corresponding OBatoms, 
    # third is the corresponding (idx) indices in the overall molecule (OBmol)
    OBatoms_in_res = []
    IDatoms_in_res = []
    NUMatoms_in_res = []
    for obatom in openbabel.OBResidueAtomIter(obres):
        OBatoms_in_res.append(obatom)
        IDatoms_in_res.append(obres.GetAtomID(obatom).strip())
        NUMatoms_in_res.append(obatom.GetIdx())

    return (IDatoms_in_res, OBatoms_in_res, NUMatoms_in_res)


def get_sidechain_atoms(obres):
    """
    TODO
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    # return a list of three lists, 
    # first list is the atoms' names, second is the corresponding OBatoms, 
    # third is the corresponding (idx) indices in the overall molecule (OBmol)
    (IDatoms_in_res, OBatoms_in_res, NUMatoms_in_res) = get_atoms(obres)
    sidechain_IDatoms = get_sidechain_atoms_names(obres)
    sidechain_OBatoms = []
    sidechain_NUMatoms = []
    for k in range(0, len(sidechain_IDatoms)):
        atomID = sidechain_IDatoms[k]
        loc = IDatoms_in_res.index(atomID)
        sidechain_OBatoms.append(OBatoms_in_res[loc])
        sidechain_NUMatoms.append(OBatoms_in_res[loc].GetIdx())

    return (sidechain_IDatoms, sidechain_OBatoms, sidechain_NUMatoms)


def get_chi_atoms(obres):
    
    """
    TODO
    
    Parameters
    ----------
    obres : OBResidue
        Query residue
    """
    
    """
    For each residuetype there is a certain set of dihedral angles in the sidechain 
    which again are defined by a certain set of atoms each.
    
    @return: returns a dictionary of arrays, describing the atoms for all the angles 
    @rtype: dictionary
    """
    residuetype = obres.GetName()
    if residuetype.upper() not in chi_atoms:
        sys.stderr.write("Warning: Unknown residuetype " + residuetype + " in residue list in MoleculeInfo.py\n")
        return
    
    return chi_atoms[residuetype.upper()]

# DEPRECATED
# def set_chi_dihedral(obmol, obres, dihedral_number, new_angle_degrees):
#     # set a side chain dihedral labeled by its number in chi_atoms dictionary, starting from 0
#     sidechain_atoms_dict = get_chi_atoms(obres)
#     
#     if (not sidechain_atoms_dict):
#         print 'residue', obres.GetName(), ' has no side chain'
#         return
#     
#     num_chi_angles = len(sidechain_atoms_dict)
# 
#     if not sidechain_atoms_dict.has_key('x' + str(dihedral_number + 1)):
#         print 'residue', obres.GetName(), ' has no chi dihedral: ', dihedral_number
#         return
#     
#     if (obres.GetName() == "PRO"):
#         print "ignoring PRO"
#         return
#     
#     sidechain_torsion_atomsID = [atomID.upper() for atomID in sidechain_atoms_dict['x' + str(dihedral_number + 1)]]
#     
#     (IDatoms_in_res, OBatoms_in_res, NUMatoms_in_res) = get_atoms(obres)
# 
#     torsion_atom_locs = []
#     torsion_obatoms = []
#     for l in range(0, 4):
#         torsion_atom_locs.append(IDatoms_in_res.index(sidechain_torsion_atomsID[l]))
#         torsion_obatoms.append(OBatoms_in_res[torsion_atom_locs[l]])
# 
#     obmol.SetTorsion(torsion_obatoms[0], torsion_obatoms[1], torsion_obatoms[2], torsion_obatoms[3], np.radians(new_angle_degrees))


if __name__ == '__main__':
    print("Starting")

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "mol2.pdb")  
    
    for j in range(0, mol.NumResidues()): 
        obres = mol.GetResidue(j)
        print(obres.GetName())
        if (getResType(obres) == 'GLY'):
            beta_atom = getBetaAtom(obres)
            print(beta_atom.GetIdx())
        
