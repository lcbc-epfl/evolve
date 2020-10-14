'''
MoleculeCreator.py

Contains a few key functions for performing mutations specifically on proteins. Extensions to more general molecules should be straightforward to implement based on the code here.

.. codeauthor:: Nicholas Browning
'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import range
from past.utils import old_div
import openbabel
from openbabel import OBResidue
from openbabel import OBAtom
from openbabel import OBMol
from openbabel import OBAtomTyper
# from openbabel import openbabel
# from openbabel.openbabel import OBResidue
# from openbabel.openbabel import OBAtom
#from openbabel.openbabel import OBMol
#from openbabel.openbabel import OBAtomTyper
import numpy as np
from . import MoleculeInfo as mi


def rotate_atoms(mol, atom_indexes, origin, rotation_axis, angle):
    """
    Applies a rotation matrix to a collection of atoms
    
    Parameters
    ----------
    mol : OBMol
        working directory
    atom_indexes : list, int
        atom indexes involved in the rotation
    origin: list
        vector used as the origin before rotation
    rotation_axis: list
        vector about which the rotation is applied
    angle: float
        angle in radians
    """
    
    rotMatrix = np.zeros(9)
    
    sn = np.sin(angle)
    cs = np.cos(angle)
    t = 1 - cs
    
    x = rotation_axis[0]
    y = rotation_axis[1]
    z = rotation_axis[2]
    
    rotMatrix[0] = t * x * x + cs
    rotMatrix[1] = t * x * y + sn * z
    rotMatrix[2] = t * x * z - sn * y
    rotMatrix[3] = t * x * y - sn * z
    rotMatrix[4] = t * y * y + cs
    rotMatrix[5] = t * y * z + sn * x
    rotMatrix[6] = t * x * z + sn * y
    rotMatrix[7] = t * y * z - sn * x
    rotMatrix[8] = t * z * z + cs
    
    for i in range (0, len(atom_indexes)):

        atom = mol.GetAtom(atom_indexes[i])
        
        xtemp = atom.GetX()
        ytemp = atom.GetY()
        ztemp = atom.GetZ()

        xtemp -= origin[0]
        ytemp -= origin[1]
        ztemp -= origin[2]

        x = xtemp * rotMatrix[0] + ytemp * rotMatrix[3] + ztemp * rotMatrix[6]
        y = xtemp * rotMatrix[1] + ytemp * rotMatrix[4] + ztemp * rotMatrix[7]
        z = xtemp * rotMatrix[2] + ytemp * rotMatrix[5] + ztemp * rotMatrix[8]

        xtemp = x
        ytemp = y
        ztemp = z
        
        xtemp += origin[0]
        ytemp += origin[1]
        ztemp += origin[2]

        atom.SetVector(xtemp, ytemp, ztemp)


def add_fragment(mol, fragment, id_start):
    """
    adds a fragment OBMol to a base mol.
    
    Parameters
    ----------
    mol : OBMol
        base mol to which fragment molecule is added
    fragment : OBMol
        fragment molecule being added
    id_start: int
        first atom in fragment will have this id after addition, may be useful for tagging/adding multiple concerted fragments.
    """
    
    prevAtoms = mol.NumAtoms()
    
    id = id_start

    for i in range (1, fragment.NumAtoms() + 1):
        frag_atom = fragment.GetAtom(i)
        frag_atom.SetId(id)
    
        new_atom = mol.NewAtom()
        new_atom.SetId(id)
        new_atom.SetVector(frag_atom.GetVector())
        new_atom.SetAtomicNum(frag_atom.GetAtomicNum())
        
        new_atom.SetType(frag_atom.GetType())
        new_atom.SetSpinMultiplicity(frag_atom.GetSpinMultiplicity())
        new_atom.SetFormalCharge(frag_atom.GetFormalCharge())
        # new_atom.SetImplicitValence(frag_atom.GetImplicitValence())
        new_atom.SetPartialCharge(frag_atom.GetPartialCharge())
        
        for data in frag_atom.GetData():
            new_atom.SetData(data)
  
        id += 1
        
    for obbond in openbabel.OBMolBondIter(fragment):
        begin_atom = obbond.GetBeginAtom()
        end_atom = obbond.GetEndAtom()
        mol.AddBond(obbond.GetBeginAtomIdx() + prevAtoms, obbond.GetEndAtomIdx() + prevAtoms, obbond.GetBondOrder(), obbond.GetFlags())
    
    
def debugAtom(obatom):
    return "AN:", obatom.GetAtomicNum(), "T:", obatom.GetType(), "idx:", obatom.GetIdx(), "id:", obatom.GetId(), "BO(1):", obatom.CountBondsOfOrder(1)


def debugBond(obbond):
    return "begin_idx:", obbond.GetBeginAtomIdx(), "end_idx:", obbond.GetEndAtomIdx(), "norm:", obbond.GetLength(), "BO:", obbond.GetBondOrder()


def getAtomByID(mol, atom_id):
    for i in range (1, mol.NumAtoms() + 1):
        if (mol.GetAtom(i).GetId() == atom_id):
            return mol.GetAtom(i)

        
def swapsidechain(settings, mol, res_index, aa_mol):
    """
    Mutates a certain residue in a protein to another amino acid. 
    
    This method first finds the CA-CB bonds for both the protein and the fragment residue. It then performs a breadth-first-search of all atoms connected to CB in the protein,
    and CA in the fragment, inclusively. These atoms are deleted, the atoms in the fragment are rotated into the correct orientation, and a new bond is made between CA of the protein and
    CB of the mutant residue. If relevant, it also updates the N-CA-CB-CG dihedral of the now mutated protein residue such that it corresponds to that of the fragment. Finally, all atom IDs
    are then modified such that they are consecutive.
    
    Parameters
    ----------
    mol : OBMol
        Base protein to perform the mutation on 
    res_index : int
        residue index to apply the mutation
    aa_mol: OBMol
        mutant residue
    """
    
    mol.BeginModify()
   
    curr = mol.GetResidue(res_index)
    
    new = aa_mol.GetResidue(0)
    print(mol, mol.NumResidues())
    print(mol.GetResidue(4))
    mol_CA = mi.getAlphaCarbon(curr)
    mol_CB = mi.getBetaAtom(curr)
    mol_N = mi.getBBNitrogen(curr)

    aa_CA = mi.getAlphaCarbon(new)
    aa_CB = mi.getBetaAtom(new)
    aa_bb_nitrogen = mi.getBBNitrogen(new)
    aa_gamma_atom = mi.getChi1DihedralAtom(new)

    frag_res = aa_CB.GetResidue()
    frag_name = frag_res.GetName()

    res_name =mi.getResType(curr)

    aa_tor = 0
    if (aa_gamma_atom is not None):
        aa_tor = aa_mol.GetTorsion(aa_gamma_atom, aa_CA, aa_CB, aa_bb_nitrogen)
    
    mol_cb_vec = np.asarray([mol_CB.GetX(), mol_CB.GetY(), mol_CB.GetZ()])
    mol_ca_vec = np.asarray([mol_CA.GetX(), mol_CA.GetY(), mol_CA.GetZ()])
    
    frag_cb_vec = np.asarray([aa_CB.GetX(), aa_CB.GetY(), aa_CB.GetZ()])
    frag_ca_vec = np.asarray([aa_CA.GetX(), aa_CA.GetY(), aa_CA.GetZ()])
    
    molbondvec = mol_cb_vec - mol_ca_vec
    fragbondvec = frag_cb_vec - frag_ca_vec
    
    mol_atoms_del = openbabel.vectorInt()
    frag_atoms_del = openbabel.vectorInt()
    
    orig_num_atoms = mol.NumAtoms()
    
    # BFS for all atoms which need to be deleted
    mol.FindChildren(mol_atoms_del, mol_CA.GetIdx(), mol_CB.GetIdx())
    aa_mol.FindChildren(frag_atoms_del, aa_CB.GetIdx(), aa_CA.GetIdx())
    
    # rotate fragment atoms such that they have the correct alignment upon adding to the protein
    if (not np.allclose(molbondvec, fragbondvec)):
        rotate_axis = np.cross(molbondvec, fragbondvec)
        rotate_axis = old_div(rotate_axis, np.linalg.norm(rotate_axis))

        dp = old_div(np.dot(molbondvec, fragbondvec), (np.linalg.norm(molbondvec) * np.linalg.norm(fragbondvec)))
        angle = np.arccos(dp)
        diffAng = -angle
        
        rotate_atoms(aa_mol, [i for i in range(1, aa_mol.NumAtoms() + 1)], frag_ca_vec, rotate_axis, diffAng)
        
        for i in range (1, aa_mol.NumAtoms() + 1):
            
            frag_atom = aa_mol.GetAtom(i)
            frag_atom.SetVector(frag_atom.GetX() - frag_ca_vec[0] + mol_ca_vec[0], frag_atom.GetY() - frag_ca_vec[1] + mol_ca_vec[1], frag_atom.GetZ() - frag_ca_vec[2] + mol_ca_vec[2])
   
    mol_atom_ids_del = [mol.GetAtom(i).GetId() for i in mol_atoms_del]
    mol_atom_ids_del.append(mol_CB.GetId())
    
    frag_atom_ids_del = [aa_mol.GetAtom(i).GetId() for i in frag_atoms_del]
    frag_atom_ids_del.append(aa_CA.GetId())
    
    # delete unnecessary atoms
    for i in range (0, len(frag_atom_ids_del)):
        
        atom = getAtomByID(aa_mol, frag_atom_ids_del[i])
        
        res = atom.GetResidue()
        
        if (res != None):
            res.RemoveAtom(atom)
        
        aa_mol.DeleteAtom(atom)
        
    for i in range (0, len(mol_atom_ids_del)):
        atom = getAtomByID(mol, mol_atom_ids_del[i])
        
        res = atom.GetResidue()
        
        if (res != None):
            res.RemoveAtom(atom)
        
        mol.DeleteAtom(atom)
        
    prev_atoms = mol.NumAtoms()
    
    # add fragment to the protein
    add_fragment(mol, aa_mol, orig_num_atoms)
    
    mol.EndModify()
    
    renum_atoms = openbabel.vectorInt(prev_atoms + aa_mol.NumAtoms())
    
    mol_CA_idx = mol_CA.GetIdx()
    
    for i in range (0, mol_CA.GetIdx()):
        renum_atoms[i] = i + 1
        
    for i in range (0, mol.NumAtoms() - prev_atoms):
        renum_atoms[mol_CA.GetIdx() + i] = prev_atoms + i + 1
        
    for i in range (0, prev_atoms - mol_CA.GetIdx()):
        renum_atoms[mol_CA.GetIdx() + mol.NumAtoms() - prev_atoms + i] = mol_CA.GetIdx() + i + 1
    
    corr_frag_cb = None
    
    for i in range(1, mol.NumAtoms() + 1):
        atom = mol.GetAtom(i)
        
        if (atom.GetId() == aa_CB.GetId()):
            corr_frag_cb = atom

    mol.AddBond(mol_CA.GetIdx(), corr_frag_cb.GetIdx(), 1)
 
    frag_res_type = mi.getResType(frag_res)

    if settings.use_res_type == True:
        curr.SetName(frag_res_type)
        pass
    else:
        curr.SetName(frag_name)
        
    for obatom in openbabel.OBResidueAtomIter(frag_res):
        frag_atom = getAtomByID(mol, obatom.GetId())
        curr.AddAtom(frag_atom)
        curr.SetAtomID(frag_atom, frag_res.GetAtomID(obatom))
    

    # For several residues in order for the amber settings to work, one needs to rename a few hydrogens which is done here.
    # if old resname = glycine we need to rename the remaining H to HA, the other has been replaced by the sidechain
    if res_name=="GLY":
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom) in ["HA2", "HA3", " HA2", " HA2 "]:
                curr.SetAtomID(obatom, "HA")
    
    if frag_res_type=="GLY" and settings.use_res_type:
        for obatom in openbabel.OBResidueAtomIter(curr):
            #print("."+curr.GetAtomID(obatom)+".")
            if curr.GetAtomID(obatom) in ["HA", " HA", " HA "] :
                curr.SetAtomID(obatom, "HA2")
    
    if frag_res_type=="HID" and settings.use_res_type:
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom)in ["HD3", "HD4", "HD5", "HD6", "HD7", "HD8"]:
                curr.SetAtomID(obatom, "HD1")

    if frag_res_type=="HID" and frag_name=="HD2" and settings.use_res_type:
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom)=="HD2":
                curr.SetAtomID(obatom, "HD1")
            if curr.GetAtomID(obatom)=="H1":
                curr.SetAtomID(obatom, "HD2")
    
    if frag_name=="HE2" and settings.use_res_type:
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom)=="HE2":
                curr.SetAtomID(obatom, "HE1")
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom)=="H1":
                curr.SetAtomID(obatom, "HE2")
    
    if frag_res_type=="HIE" and frag_name!="HE2" and settings.use_res_type:
        for obatom in openbabel.OBResidueAtomIter(curr):
            if curr.GetAtomID(obatom) in ["HE3", "HE4", "HE5", "HE6", "HE7", "HE8"]:
                curr.SetAtomID(obatom, "HE1")
               




    alpha_carbon = mi.getAlphaCarbon(curr)
    bb_nitrogen = mi.getBBNitrogen(curr)
    beta_atom = mi.getBetaAtom(curr)
    chi_atom = mi.getChi1DihedralAtom(curr)
    
    if (chi_atom is not None and aa_gamma_atom is not None):
        mol.SetTorsion(bb_nitrogen, alpha_carbon, beta_atom, chi_atom, aa_tor * (old_div(np.pi, 180.0)))
  
    # need to renumber IDs now to be consecutive
    mol.RenumberAtoms(renum_atoms)
    
    for i in range (1, mol.NumAtoms() + 1) :
        mol.GetAtom(i).SetId(i - 1)


        
    
if __name__ == '__main__':
    pass
    
