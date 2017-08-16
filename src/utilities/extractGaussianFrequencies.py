'''
main.py

@author: Nicholas Browning
'''
import numpy as np

def getFrequencyDict(gaussian_out_file, num_atoms):
    '''atomic_numbers, frequency_index -> {frequency, delta_xzys}'''
    #                     1                      2                      3
    #                      A                      A                      A
    # Frequencies --     10.4859                12.2013                14.7020
    # Red. masses --      5.4037                 5.3983                 5.5738
    # Frc consts  --      0.0004                 0.0005                 0.0007
    # IR Inten    --      3.4759                 0.3909                 0.6587
    #  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
    
    frequency_dict = {}
    
    file_p = open(gaussian_out_file, 'r')

    file_lines = file_p.readlines()

    file_p.close()

    i = 0

    while i < len(file_lines):
        line = file_lines[i]
    
        if "Frequencies --" in line:
            # appears in packs of 3
            frequencies = [float(v) for v in line.split()[2:]] 
            freq_indexes = [int(v) for v in file_lines[i - 2].split()]
    
            for k in xrange(0, 3):
                limits = [[2, 5], [5, 8], [8, 11]]
                atomic_numbers = []
                delta_positions = []
            
                for j in xrange(i + 5, i + 5 + num_atoms):
                    curr_line = file_lines[j].split()
                    freq_xyz = curr_line[limits[k][0]:limits[k][1]]
                    atomic_numbers.append(int(curr_line[1]))
                    delta_positions.append([float(v) for v in freq_xyz])
                
                frequency_dict[freq_indexes[k]] = (np.asarray(frequencies[k]), np.asarray(delta_positions))
    
            i = j
        i += 1
    return np.asarray(atomic_numbers, dtype=int), frequency_dict


def getVibratingAtoms(atom_numbers, delta_positions):
    norms = getVibrationNorms(delta_positions)
    indexes = np.argsort(norms)[::-1]
    non_zero_norm_indexes = np.where(norms[indexes] > 0)
    non_zero_norm_indexes = non_zero_norm_indexes[0]

    # for i in non_zero_norm_indexes:
        # print i, indexes[i], atom_numbers[indexes][i], norms[indexes][i], delta_positions[indexes][i]
    return indexes[non_zero_norm_indexes]
        
def getVibrationNorms(delta_positions):
    return np.linalg.norm(delta_positions, axis=1)
    

    
