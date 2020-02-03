'''
constants file

.. codeauthor:: Nicholas Browning
.. codeauthor:: Simon Duerr dev@simonduerr.eu

'''
import numpy as np

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
charge2element = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9:'F', 32:'S'}

bohr2angstrom = 0.5291772083
angstrom2bohr = 1.88972613386246

rotamers = ['G01', 'A01', 'R01', 'R02', 'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', \
            'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18', 'R19', \
            'R20', 'R21', 'R22', 'R23', 'R24', 'R25', 'R26', 'R27', 'R28', 'D01', \
            'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'N01', \
            'N02', 'N03', 'N04', 'N05', 'N06', 'N07', 'C01', 'C02', 'C03', 'E01', \
            'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09', 'E10', 'E11', \
            'E12', 'E13', 'E14', 'Q01', 'Q02', 'Q03', 'Q04', 'Q05', 'HD1', \
            'HD2', 'HD3', 'HD4', 'HD5', 'HD6', 'HD7', 'HD8', 'HE1', 'HE2', 'HE3', \
            'HE4', 'HE5', 'HE6', 'HE7', 'HE8', 'HP1', 'HP2', 'HP3', 'HP4', 'HP5', \
            'HP6', 'HP7', 'HP8', 'I01', 'I02', 'I03', 'I04', 'I05', 'L01', 'L02', \
            'L03', 'L04', 'K01', 'K02', 'K03', 'K04', 'K05', 'K06', 'K07', 'K08', \
            'K09', 'K10', 'K11', 'K12', 'K13', 'K14', 'K15', 'K16', 'K17', 'K18', \
            'K19', 'K20', 'K21', 'K22', 'K23', 'K24', 'K25', 'K26', 'K27', 'K28', \
            'K29', 'K30', 'K31', 'K32', 'K33', 'K34', 'K35', 'K36', 'K37', 'K38', \
            'M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', \
            'M11', 'M12', 'M13', 'P01', 'P02', 'P03', 'P04', 'S01', 'S02', 'S03', \
            'T01', 'T02', 'T03', 'W01', 'W02', 'W03', 'W04', 'W05', 'W06', 'W07', \
            'Y01', 'Y02', 'Y03', 'Y04', 'V01', 'V02', 'V03'
            ]

selected_rotamers = rotamers

# ARG, GLU, TRP, HID, HIP
rotamer_types = [ 'GLY', 'ALA', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', \
            'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', \
            'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ASP', \
            'ASP', 'ASP', 'ASP', 'ASP', 'ASH', 'ASH', 'ASH', 'ASH', 'ASH', 'ASN', \
            'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'CYS', 'CYS', 'CYS', 'GLU', \
            'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLH', 'GLH', 'GLH', 'GLH', \
            'GLH', 'GLH', 'GLH', 'GLN', 'GLN', 'GLN', 'GLN', 'GLN', 'HID', \
            'HID', 'HID', 'HID', 'HID', 'HID', 'HID', 'HID', 'HIE', 'HIE', 'HIE', \
            'HIE', 'HIE', 'HIE', 'HIE', 'HIE', 'HIP', 'HIP', 'HIP', 'HIP', 'HIP', \
            'HIP', 'HIP', 'HIP', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'LEU', 'LEU', \
            'LEU', 'LEU', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', \
            'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', 'LYS', \
            'LYS', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', \
            'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', 'LYN', \
            'MET', 'MET', 'MET', 'MET', 'MET', 'MET', 'MET', 'MET', 'MET', 'MET', \
            'MET', 'MET', 'MET', 'PHE', 'PHE', 'PHE', 'PHE', 'SER', 'SER', 'SER', \
            'THR', 'THR', 'THR', 'TRP', 'TRP', 'TRP', 'TRP', 'TRP', 'TRP', 'TRP', \
            'TYR', 'TYR', 'TYR', 'TYR', 'VAL', 'VAL', 'VAL'
            ]
selected_rotamer_types = rotamer_types

allowed_residue_types = ['GLY', 'ALA', 'ARG', 'ASP', 'ASH', 'ASN', 'CYS', 'GLU', 'GLH', 'GLN', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'LYN', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# 10, 50, 80
energies = {'ALA': (13.03, 12.14 , 12.06), 'ARG': (-170.77, -176.22 , -176.73), 'ASH': (-43.74, -45.64, -45.82), 'ASP': (-52.27, -58.92, -59.55), 'ASN': (-72.72, -74.29, -74.44), \
            'CYS': (13.55, 12.58, 12.49), 'GLH': (-36.86, -38.93, -39.13), 'GLU': (-47.00, -53.61, -54.30), 'GLN': (-56.93, -58.79, -58.97), 'GLY':(7.11, 6.15, 6.06), \
            'HID': (12.43, 10.83, 10.68), 'HIE':(8.74, 7.22, 7.08), 'HIP': (22.66, 17.17, 16.62), 'ILE': (10.74, 9.93, 9.85), 'LEU': (-10.04, -10.88, -10.96), \
            'LYN': (-4.36, -5.73, -5.85), 'LYS': (-9.60, -15.72, -16.29), 'MET': (8.36, 7.47, 7.29), 'PHE': (11.05, 10.04, 9.95), 'SER': (-0.71, -2.32, -2.47), \
            'THR': (-19.37, -20.45, -20.59), 'TRP':(14.43, 13.07, 12.94), 'TYR': (-13.28, -14.89, -15.05), 'VAL':(-7.58, -8.41, -8.49) }

def subselect_rotamers(type_list):
    indices = [i for i, x in enumerate(rotamer_types) if (x in type_list)]
    return [rotamers[i] for i in indices]

def subselect_rotamer_types(type_list):
    indices = [i for i, x in enumerate(rotamer_types) if (x in type_list)]
    return [rotamer_types[i] for i in indices]

def subselect_rotamer_energy_dict(type_list):
    new_dict = {}
    
    for i in type_list:
        new_dict[i] = energies[i]
        
    return new_dict
