'''
main.py

@author: Nicholas Browning
'''

import openbabel
from openbabel import OBResidue
from openbabel import OBAtom
from openbabel import OBMol
from openbabel import OBAtomTyper
import numpy as np
import MoleculeInfo as mi


def rotate_atoms(mol, atom_indexes, origin, rotation_axis, angle):
    
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
    
    for i in xrange (0, len(atom_indexes)):

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
    
    prevAtoms = mol.NumAtoms()
    
    id = id_start
    
    # copy atomic data
    for i in xrange (1, fragment.NumAtoms() + 1):
        frag_atom = fragment.GetAtom(i)
        frag_atom.SetId(id)
    
        new_atom = mol.NewAtom()
        new_atom.SetId(id)
        new_atom.SetVector(frag_atom.GetVector())
        new_atom.SetAtomicNum(frag_atom.GetAtomicNum())
        
        new_atom.SetType(frag_atom.GetType())
        new_atom.SetSpinMultiplicity(frag_atom.GetSpinMultiplicity())
        new_atom.SetFormalCharge(frag_atom.GetFormalCharge())
        new_atom.SetImplicitValence(frag_atom.GetImplicitValence())
        new_atom.SetPartialCharge(frag_atom.GetPartialCharge())
    
        
        
        for data in frag_atom.GetData():
            new_atom.SetData(data)
            
        # print "Added atom:", debugAtom(new_atom)
        
        id += 1
        
    # copy bond data    
    for obbond in openbabel.OBMolBondIter(fragment):
        begin_atom = obbond.GetBeginAtom()
        end_atom = obbond.GetEndAtom()
        
        mol.AddBond(obbond.GetBeginAtomIdx() + prevAtoms, obbond.GetEndAtomIdx() + prevAtoms, obbond.GetBO(), obbond.GetFlags())
    
    
def debugAtom(obatom):
    return "AN:", obatom.GetAtomicNum(), "T:", obatom.GetType(), "idx:", obatom.GetIdx(), "id:", obatom.GetId(), "BO(1):", obatom.CountBondsOfOrder(1)

def debugBond(obbond):
    return "begin_idx:", obbond.GetBeginAtomIdx(), "end_idx:", obbond.GetEndAtomIdx(), "norm:", obbond.GetLength(), "BO:", obbond.GetBO()


def getAtomByID(mol, atom_id):
    for i in xrange (1, mol.NumAtoms() + 1):
        if (mol.GetAtom(i).GetId() == atom_id):
            return mol.GetAtom(i)
        
def swapsidechain (mol, res_index, aa_mol):
    
    mol.BeginModify()
    
    curr = mol.GetResidue(res_index)
    
    # print curr.GetName()
    
    new = aa_mol.GetResidue(0)
    
    mol_CA = mi.getAlphaCarbon(curr)
    mol_CB = mi.getBetaAtom(curr)
    mol_N = mi.getBBNitrogen(curr)

    # print "MOL_CA", debugAtom(mol_CA)
    # print "MOL_CB", debugAtom(mol_CB)
    
    aa_CA = mi.getAlphaCarbon(new)
    aa_CB = mi.getBetaAtom(new)
    aa_bb_nitrogen = mi.getBBNitrogen(new)
    aa_gamma_atom = mi.getChi1DihedralAtom(new)
    
    
    frag_res = aa_CB.GetResidue()
    frag_name = frag_res.GetName()
  
    aa_tor = aa_mol.GetTorsion(aa_gamma_atom, aa_CA, aa_CB, aa_bb_nitrogen)
    
    # print aa_tor
    
    # print "FRAG_CA", debugAtom(aa_CA)
    # print "FRAG_CB", debugAtom(aa_CB)
    
    
    mol_cb_vec = np.asarray([mol_CB.GetX(), mol_CB.GetY(), mol_CB.GetZ()])
    mol_ca_vec = np.asarray([mol_CA.GetX(), mol_CA.GetY(), mol_CA.GetZ()])
    
    frag_cb_vec = np.asarray([aa_CB.GetX(), aa_CB.GetY(), aa_CB.GetZ()])
    frag_ca_vec = np.asarray([aa_CA.GetX(), aa_CA.GetY(), aa_CA.GetZ()])
    
    molbondvec = mol_cb_vec - mol_ca_vec
    fragbondvec = frag_cb_vec - frag_ca_vec
    
    mol_atoms_del = openbabel.vectorInt()
    frag_atoms_del = openbabel.vectorInt()
    
    orig_num_atoms = mol.NumAtoms()
    # print orig_num_atoms, aa_mol.NumAtoms()
    
    mol.FindChildren(mol_atoms_del, mol_CA.GetIdx(), mol_CB.GetIdx())
    aa_mol.FindChildren(frag_atoms_del, aa_CB.GetIdx(), aa_CA.GetIdx())
    
    if (not np.allclose(molbondvec, fragbondvec)):
        rotate_axis = np.cross(molbondvec, fragbondvec)
        rotate_axis = rotate_axis / np.linalg.norm(rotate_axis)

        dp = np.dot(molbondvec, fragbondvec) / (np.linalg.norm(molbondvec) * np.linalg.norm(fragbondvec))
        angle = np.arccos(dp)
        diffAng = -angle
        
        rotate_atoms(aa_mol, [i for i in xrange(1, aa_mol.NumAtoms() + 1)], frag_ca_vec, rotate_axis, diffAng)
        
        for i in xrange (1, aa_mol.NumAtoms() + 1):
            
            frag_atom = aa_mol.GetAtom(i)
            frag_atom.SetVector(frag_atom.GetX() - frag_ca_vec[0] + mol_ca_vec[0], frag_atom.GetY() - frag_ca_vec[1] + mol_ca_vec[1], frag_atom.GetZ() - frag_ca_vec[2] + mol_ca_vec[2])
   
    mol_atom_ids_del = [mol.GetAtom(i).GetId() for i in mol_atoms_del]
    mol_atom_ids_del.append(mol_CB.GetId())
    
    
    frag_atom_ids_del = [aa_mol.GetAtom(i).GetId() for i in frag_atoms_del]
    frag_atom_ids_del.append(aa_CA.GetId())
    
    
    for i in xrange (0, len(frag_atom_ids_del)):
        
        atom = getAtomByID(aa_mol, frag_atom_ids_del[i])
        
        res = atom.GetResidue()
        
        if (res != None):
            res.RemoveAtom(atom)
        
        aa_mol.DeleteAtom(atom)
        
    for i in xrange (0, len(mol_atom_ids_del)):
        atom = getAtomByID(mol, mol_atom_ids_del[i])
        
        res = atom.GetResidue()
        
        if (res != None):
            res.RemoveAtom(atom)
        
        mol.DeleteAtom(atom)
        
    prev_atoms = mol.NumAtoms()
    

    # print mol.NumAtoms(), aa_mol.NumAtoms()
    
    # print "DEBUGGING MOL BONDS:"
    # for i in xrange (0, mol.NumBonds()):
    #    print debugBond(mol.GetBond(i))
        
    # print "DEBUGGING FRAG BONDS:"
    # for i in xrange (0, aa_mol.NumBonds()):
    #    print debugBond(aa_mol.GetBond(i))
    
    # print "Adding Fragment"
    add_fragment(mol, aa_mol, orig_num_atoms)
    mol.EndModify()
    renum_atoms = openbabel.vectorInt(prev_atoms + aa_mol.NumAtoms())
    
    mol_CA_idx = mol_CA.GetIdx()
    
    for i in xrange (0, mol_CA.GetIdx()):
        renum_atoms[i] = i + 1
        
    for i in xrange (0, mol.NumAtoms() - prev_atoms):
        renum_atoms[mol_CA.GetIdx() + i] = prev_atoms + i + 1
        
    for i in xrange (0, prev_atoms - mol_CA.GetIdx()):
        renum_atoms[mol_CA.GetIdx() + mol.NumAtoms() - prev_atoms + i] = mol_CA.GetIdx() + i + 1
    
    # print mol.NumAtoms()
    
    # print "DEBUGGING MOL BONDS:"
    # for i in xrange (0, mol.NumBonds()):
    #    print debugBond(mol.GetBond(i))
    
    # print "DEBUG MOL ATOMS:"
    # for i in xrange(1, mol.NumAtoms() + 1):
     #   atom = mol.GetAtom(i)
     #   print debugAtom(atom)
        
    # print mol_CA, mol_CA.GetId() 
    # print aa_CB, aa_CB.GetId()
    
    corr_frag_cb = None
    
    for i in xrange(1, mol.NumAtoms() + 1):
        atom = mol.GetAtom(i)
        
        if (atom.GetId() == aa_CB.GetId()):
            corr_frag_cb = atom
                
    # print corr_frag_cb
    
    # print "corresponding frag cb:", debugAtom(corr_frag_cb)
    # print "mol ca: ", debugAtom(mol_CA)
    # print "Adding new bond"
    

    mol.AddBond(mol_CA.GetIdx(), corr_frag_cb.GetIdx(), 1)
    
    
    
    # print "DEBUGGING MOL BONDS:"
    # for i in xrange (0, mol.NumBonds()):
    #    print debugBond(mol.GetBond(i))
    
    # print "Debug Frag Res:"
    # for obatom in openbabel.OBResidueAtomIter(frag_res):
        # print debugAtom(obatom)
        
     
    # print "DEBUG FRAG ATOMS:"
    # for i in xrange(1, aa_mol.NumAtoms() + 1):
     #   atom = aa_mol.GetAtom(i)
        # print debugAtom(atom)
        
           
    # print "DEBUG MOL ATOMS:"
    # for i in xrange(1, mol.NumAtoms() + 1):
    #    atom = mol.GetAtom(i)
        # print debugAtom(atom)
        
    # for data in frag_res.GetData():
     #   print data
    
    
    frag_res_type = mi.getResType(frag_res)
    
    curr.SetName(frag_res_type)

    # print frag_res.GetNumAtoms()
        
    for obatom in openbabel.OBResidueAtomIter(frag_res):
        # print obatom.GetId(), obatom.GetType()
        frag_atom = getAtomByID(mol, obatom.GetId())
        # print frag_atom.GetId(), frag_atom.GetType()
        curr.AddAtom(frag_atom)
        # print curr.GetAtomID(frag_atom), frag_res.GetAtomID(obatom)
        curr.SetAtomID(frag_atom, frag_res.GetAtomID(obatom))
        # print curr.GetAtomID(frag_atom), frag_res.GetAtomID(obatom)
        
    # print curr.GetNumAtoms()
        
    alpha_carbon = mi.getAlphaCarbon(curr)
    bb_nitrogen = mi.getBBNitrogen(curr)
    beta_atom = mi.getBetaAtom(curr)
    chi_atom = mi.getChi1DihedralAtom(curr)
    
    
    
    # print "tor", alpha_carbon, bb_nitrogen, beta_atom, chi_atom
    # print debugAtom(alpha_carbon)
    # print debugAtom(bb_nitrogen)
    # print debugAtom(beta_atom)
    # print debugAtom(chi_atom)
    
    # print "torsion update"
    mol.SetTorsion(bb_nitrogen, alpha_carbon, beta_atom, chi_atom, aa_tor * (np.pi / 180.0))

    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdb")
    obConversion.WriteFile(mol, "test2.pdb")
    # print "Renumbering Atoms"
    # print "renumber"
    mol.RenumberAtoms(renum_atoms)
    # print "renumbered"
    for i in xrange (1, mol.NumAtoms() + 1) :
        mol.GetAtom(i).SetId(i - 1)
        
    
    
if __name__ == '__main__':
    print "Starting"

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    
    obtype = openbabel.OBAtomTyper()
    obConversion.ReadFile(mol, "gpgg.pdb") 
    
        
    frag = openbabel.OBMol()
    obConversion.SetInAndOutFormats("mol2", "pdb")
    obConversion.ReadFile(frag, "/lcbcdata/nbrownin/Rotamers_mol2/ASH/D07.mol2") 
    
    frag2 = openbabel.OBMol()
    obConversion.SetInAndOutFormats("mol2", "pdb")
    obConversion.ReadFile(frag2, "/lcbcdata/nbrownin/Rotamers_mol2/HID/HD3.mol2") 
    

    print "first swap"
    
    swapsidechain(mol, 1, frag)
    
    print "second swap"
    swapsidechain(mol, 3, frag2)
    print "WRITING OUT"
    obConversion.WriteFile(mol, "test2.pdb")
    
    