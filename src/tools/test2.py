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

    directory = '/share/lcbcsrv5/lcbcdata/nbrownin/EVOLVE_A20/SINGLE_ROT_HELIX/'
    folders = ['e10', 'e50', 'e80']
    dielectric = [0, 1, 2]
    labels = ['10', '50', '80']
 
    def get_stabilities(file_name, dielectric):
        fp1 = open(file_name, 'r')
        file_lines = fp1.readlines()
        fp1.close()
    
        helix_energies_str = file_lines[2:len(file_lines): 2]
        helix_names_str = file_lines[1:len(file_lines): 2]
        for i in xrange(len(helix_names_str)):
            helix_names_str[i] = helix_names_str[i].split(':')[1].replace('.pdb', '').replace(' ', '').replace('\n', '')
        
        helix_energies = np.zeros(len(helix_energies_str))
        for i, v in enumerate(helix_energies_str):
            helix_energies[i] = np.float(v)
        
        stability = np.zeros(len(cnts.rotamers))
        
        
        res_dict_min_stability = {}
        res_dict_min_rot = {}
        for i in xrange(0, len(cnts.rotamer_types)):
            if (not res_dict_min_stability.has_key(cnts.rotamer_types[i])):
                res_dict_min_stability[cnts.rotamer_types[i]] = 999.0
                res_dict_min_rot[cnts.rotamer_types[i]] = None
            
        
        for i, v in enumerate(helix_names_str):
            
            rotamer_energy = cnts.energies[cnts.rotamer_types[cnts.rotamers.index(v)]][dielectric]
            
            stability[i] = helix_energies[i] - helix_energies[0] + ((8 * cnts.energies['ALA'][dielectric]) - (8 * rotamer_energy))
            # print helix_names_str[i], stability[i], cnts.rotamers.index(v), rotamer_energy
            
        for i, v in enumerate(helix_names_str):
            rot_type = cnts.rotamer_types[cnts.rotamers.index(v)]
            if (res_dict_min_stability[rot_type] == 999.0 or stability[i] < res_dict_min_stability[rot_type]):
                res_dict_min_stability[rot_type] = stability[i]
                res_dict_min_rot[rot_type] = v
                     
        return res_dict_min_stability, res_dict_min_rot
    

    import matplotlib.pyplot as plt
    sns.set(font_scale=2.0, rc={'text.usetex' : True})
    sns.set_style("ticks")
    sns.set_palette('Set2')

    fig, ax = plt.subplots(figsize=(12, 6))
    bw = 0.32
    
    rv = 0
    for i, v in enumerate(folders):
        res_dict_min_stability, res_dict_min_rot = get_stabilities(directory + v + '/homo_helix.dat', dielectric[i])
        
        
        print res_dict_min_stability.keys()
        x = np.arange(len(res_dict_min_stability.keys()))
        
        set_xticks(ax)
        print res_dict_min_stability
        ax.bar(x + rv, res_dict_min_stability.values(), width=bw, label=r'$\epsilon=' + labels[i] + '$')
        ax.set_ylim(-50.0, 150.0)
        for j, height in enumerate(res_dict_min_stability.values()):
            
            chosen_va = 'bottom' 
            offset = 2
            if (height < 0):
                chosen_va = 'top'
                offset = -2
                
            ax.text(x[j] + rv, height + offset, res_dict_min_rot.values()[j],
            ha='center', va=chosen_va, fontsize=10, rotation=90)
            
        rv += bw
        
    ax.set_xticks(x + bw)
    ax.set_xticklabels(res_dict_min_stability.keys(), rotation=90)
    
    
    ax.set_ylabel(r'fitness $(\text{kcal mol}^{-1})$')
    plt.legend()
    
    plt.tight_layout()
    #plt.show()
    plt.savefig('rotamer_hist_all_die.png', dpi=300)
    


