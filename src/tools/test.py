'''
main.py

@author: Nicholas Browning
'''

import sys
import argparse 
import openbabel
from src import outprocesses as op
import subprocess
from src import MoleculeInfo as mi
from src import MoleculeCreator as mc
from src import constants as cnts
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def set_xticks(ax):
    ax.tick_params(which='major', direction='in', length=7, width=1.5)

params = {
          'text.latex.preamble': [r"\usepackage{upgreek}",
                                  r"\usepackage[nice]{units}", r"\usepackage{amsmath}"]}
plt.rcParams.update(params)

def helical_stability(initial_mol, mol, mutated_indexes, DIE):
    from src import constants
    from src import MoleculeInfo as mi
     
    add = 0.0
    negate = 0.0
        
    for j in xrange (0, len(mutated_indexes)):
        res = mi.getResType(mol.GetResidue(int(mutated_indexes[j])))
        add += constants.energies['ALA'][DIE]
        negate += constants.energies[res][DIE]
   
    return (add - negate)


if __name__ == '__main__':
    
    # diff = helical_stability(initial_mol, mol, args.mutate_indexes, args.dielectric)
    
    # obConversion.WriteFile(mol, args.out_molecule)
    
    # print opt_e - initial_e + diff
    
    vacuum_lib = {
                  'ASH':-23.164,
                  'ASN':-55.627,
                  'CYS': 24.136,
                  'GLH':-16.046,
                  'GLN':-38.756,
                  'ILE': 19.622,
                  'LEU':-0.74,
                  'LYN': 10.42,
                  'MET': 18.15,
                  'PHE': 22.104,
                  'SER': 16.158,
                  'THR':-3.08,
                  'TRP': 28.982,
                  'TYR': 4.69,
                  'VAL': 1.574}
    
    targets = [ 'G01', 'A01', 'D06', 'N05', 'C01', 'E08', 'Q03', 'I01', 'L01', 'K20', 'M06', 'P02', 'S03', 'T01', 'W02', 'Y02', 'V01']
    
    allowed_types = ['GLY', 'ALA', 'ASH', 'ASN', 'CYS', 'GLH', 'GLN', 'ILE', 'LEU', 'LYN', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    
    rotamers = cnts.subselect_rotamers(allowed_types)
    rotamer_types = cnts.subselect_rotamer_types(allowed_types)
    
    fp1 = open('homo_helix.dat', 'r')
    file_lines = fp1.readlines()
    fp1.close()

    helix_energies_str = file_lines[2:len(file_lines): 2]
    helix_names_str = file_lines[1:len(file_lines): 2]
    for i in xrange(len(helix_names_str)):
        helix_names_str[i] = helix_names_str[i].split(':')[1].replace('.pdb', '').replace(' ', '').replace('\n', '')
        
    fp1 = open('rot_ener_vacuum.dat', 'r')
    file_lines = fp1.readlines()
    fp1.close()

    rotamer_energies_str = file_lines[1:len(file_lines): 2]
    rotamer_names_str = file_lines[0:len(file_lines): 2]

    for i in xrange(len(rotamer_names_str)):
        rotamer_names_str[i] = rotamer_names_str[i].split(' ')[1].replace('.pdb', '').replace(' ', '').replace('\n', '')
    
    print helix_names_str
    print rotamer_names_str
    
    rotamer_energies = np.zeros(len(rotamer_energies_str))
    for i, v in enumerate(rotamer_energies_str):
        rotamer_energies[i] = np.float(v)
        
    helix_energies = np.zeros(len(helix_energies_str))
    for i, v in enumerate(helix_energies_str):
        helix_energies[i] = np.float(v)
    
    print helix_energies
    print helix_names_str
    stability = np.zeros(len(targets))
    
    dft_stabilities = np.asarray([14.81, 0.00, 5.63, 26.26, 22.45, -40.61, -35.63 , -21.20, -27.94, -3.816, -28.87, -25.29, 29.75, 11.96, -27.41, -24.45, -7.72])
    
    for i, v in enumerate(targets):
        target_idx = 0
        for j in xrange(0, len(helix_names_str)):
            if (helix_names_str[j] == v):
                target_idx = j
                break

        rotamer_energy = rotamer_energies[rotamer_names_str.index(helix_names_str[target_idx])]
        
        stability[i] = helix_energies[target_idx] - helix_energies[0] + ((8 * rotamer_energies[0]) - (8 * rotamer_energy))
        print helix_names_str[target_idx], stability[i]
    
    import matplotlib.pyplot as plt
    sns.set(font_scale=2.0, rc={'text.usetex' : True})
    sns.set_style("ticks")
    sns.set_palette('Set2')

    print dft_stabilities
    print stability
    bar_width = 0.35
    x = np.arange(len(targets))
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    set_xticks(ax)
    ax.bar(x, dft_stabilities, width=bar_width, alpha=0.5, label='DFT')
    ax.bar(x + bar_width, stability, width=bar_width, alpha=0.5, label='MM')
    # sns.barplot(targets, )
    ax.set_xticks(x + (bar_width / 2))
    ax.set_xticklabels(allowed_types, rotation=90)
    ax.set_ylabel(r'fitness $\text{ (kcal mol)}^{-1}$')
    # plt.show()
    plt.legend()
    plt.savefig('best_MM.png', dpi=300)
