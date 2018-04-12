'''
main.py

@author: Nicholas Browning
'''

import sys

def parseMMBPSA_total_decomp(mmpbsa_file_path, num_res):
    mmpbsa_file = open(mmpbsa_file_path, 'r')
    
    lines = mmpbsa_file.readlines()
    
    mmpbsa_file.close()
    
    start = -1
    
    for i in range(len(lines)):
        
        if ("Total Energy Decomposition" in lines[i]):
            start = i
            break
    
    
    curr = i + 3
    end = i + 3 + (num_res)
    print curr, end
    
    
    # res name = 0
    # internal_energy = 1
    # VdW = 4
    # electrostatics = 7
    # polar solvation = 10
    # non-polar solvation = 13
    # total energy = 16
    
    
    count = 0
    return_dict = {}
    while (curr < end):
        data = lines[curr].split(',')
        print (data[0], data[1], data[4], data[7], data[10], data[13], data[16])
        
        return_dict[count] = (data[1], data[4], data[7], data[10], data[13], data[16])

        curr += 1
        count += 1
    
    return return_dict

if __name__ == "__main__":
    values_dict = parseMMBPSA_total_decomp(sys.argv[1], 20)
    
    print values_dict[0][5]
