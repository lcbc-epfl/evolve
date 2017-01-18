'''

optimise.py

@author: Nicholas Browning
'''

import numpy as np
import matplotlib.pyplot as plt

'''
general, gly, prepro, pro
'''

import os
import sys

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def read_input_file(inputFilePath):
    import ConfigParser
    settings = ConfigParser.ConfigParser()
    settings.read(inputFilePath)
    return settings

def loadold():
    general_data = np.fromfile('./psiphi/rama500-general-nosec.data', sep=' ')
    gly_data = np.fromfile('./psiphi/rama500-gly-sym-nosec.data', sep=' ')
    prepro_data = np.fromfile('./psiphi/rama500-prepro.data', sep=' ')
    pro_data = np.fromfile('./psiphi/rama500-pro.data', sep=' ')

    # x = data[:, 0]
    # y = data[:, 1]
    
    general_data = general_data.reshape(len(general_data) / 3, 3)

    general_probability = general_data[:, 2]
  
    # general_probability = general_probability / max(general_probability)

    N = int(len(general_probability) ** .5)
    general_probability = general_probability.reshape(N, N, order='C')

    gly_data = gly_data.reshape(len(gly_data) / 3, 3)
    gly_probability = gly_data[:, 2]
    # gly_probability = gly_probability / max(gly_probability)
    gly_probability = gly_probability.reshape(N, N, order='C')
    
    prepro_data = prepro_data.reshape(len(prepro_data) / 3, 3)
    prepro_probability = prepro_data[:, 2]
    # prepro_probability = prepro_probability / max(prepro_probability)
    prepro_probability = prepro_probability.reshape(N, N, order='C')
    
    pro_data = pro_data.reshape(len(pro_data) / 3, 3)
    pro_probability = pro_data[:, 2]
    # pro_probability = pro_probability / max(pro_probability)
    pro_probability = pro_probability.reshape(N, N, order='C')
    
    return [general_probability, gly_probability, prepro_probability, pro_probability]

def loadDefaultProbabilityDistributions(distribution_paths='/psiphi/rama500-general-nosec.data /psiphi/rama500-gly-sym-nosec.data /psiphi/rama500-prepro.data /psiphi/rama500-pro.data'):
    script_path = get_script_path()
    distributions = []
    for distribution_data in distribution_paths.split():
        data = np.fromfile(script_path + distribution_data, sep=' ')
        data = data.reshape(len(data) / 3, 3)
        data_probability = data[:, 2]
        N = int(len(data_probability) ** .5)
        data_probability = data_probability.reshape(N, N, order='F')
        distributions.append(data_probability)
    return distributions
    

def plotProbabilityDistribution(probabilityDistribution):
    plt.imshow(probabilityDistribution, cmap=plt.get_cmap("rainbow"), extent=[-180, 180, -180, 180], origin="lower")
    plt.colorbar()
    plt.show()
    
def plotProbabilityDistributionAndSolution(probabilityDistributions, probabilityDistributionPointers, solution):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    axes = [ax1, ax2, ax3, ax4]
    
    ax1.imshow(probabilityDistributions[0], cmap=plt.get_cmap("rainbow"), origin="lower")
    ax1.set_title('General')
    ax1.set_xlabel('$\Phi$')
    ax1.set_ylabel('$\Psi$')
    ax1.set_xlim([0, 180])
    ax1.set_ylim([0, 180])
    
    ax2.imshow(probabilityDistributions[1], cmap=plt.get_cmap("rainbow"), origin="lower")
    ax2.set_title('GLY')
    ax2.set_xlabel('$\Phi$')
    ax2.set_ylabel('$\Psi$')
    ax2.set_xlim([0, 180])
    ax2.set_ylim([0, 180])
    
    ax3.imshow(probabilityDistributions[2], cmap=plt.get_cmap("rainbow"), origin="lower")
    ax3.set_title('PREPRO')
    ax3.set_xlabel('$\Phi$')
    ax3.set_ylabel('$\Psi$')
    ax3.set_xlim([0, 180])
    ax3.set_ylim([0, 180])
    
    ax4.imshow(probabilityDistributions[3], cmap=plt.get_cmap("rainbow"), origin="lower")
    ax4.set_title('PRO')
    ax4.set_xlim([0, 180])
    ax4.set_ylim([0, 180])
    ax4.set_xlabel('$\Phi$')
    ax4.set_ylabel('$\Psi$')
    
    for i, sol in enumerate(solution):
        axes[probabilityDistributionPointers[i]].scatter((180 + sol[0]) / 2, (180 + sol[1]) / 2, c='black', s=15)
    plt.tight_layout()
    plt.show()
    
def doTest():
    pass

def doMain():
    import argparse
    from distutils.util import strtobool

    parser = argparse.ArgumentParser()
    parser.add_argument("-input_file", help="input file for montecarlo dihedral optimisation", type=str)
    args = parser.parse_args()
    
    settings = read_input_file(args.input_file)
    print  settings
    residue_pointers = [int (v) for v in  settings.get('SETTINGS', 'RESIDUE_POINTERS').split()]
    plot_optimised_distribution = strtobool(settings.get('SETTINGS', 'PLOT_OPTIMISED_DISTRIBUTION'))
    
    
    distribution_path = settings.get('DISTRIBUTION', 'DISTRIBUTION_FILES')
        
    probs = loadDefaultProbabilityDistributions(distribution_paths=distribution_path)
    
    num_phipsi = len(residue_pointers)
    
    initial_solution = np.random.uniform(low=-180, high=180, size=(num_phipsi, 2))
    
    print "Initial Solution:"
    print initial_solution
    
    move_lows = [-180] * num_phipsi
    move_highs = [180] * num_phipsi
    
    final_solution = optimisePsiPhi(probs, residue_pointers, initial_solution, int(settings.get('SETTINGS', 'NUM_MC_STEPS')), move_lows, move_highs)
    
    print "Final Solution:"
    print final_solution
    
    
    if (plot_optimised_distribution):
        plotProbabilityDistributionAndSolution(probs, residue_pointers, final_solution)
    
    if settings.has_section('MOLBUILD'):
        from api import RotamerCreator as rc
        molbuild_residue_numbers = settings.get('MOLBUILD', 'MOLBUILD_RESIDUE_NUMBERS').split()
        molbuild_input = settings.get('MOLBUILD', 'MOLBUILD_INPUT')
        molbuild_output = settings.get('MOLBUILD', 'MOLBUILD_OUTPUT')
        # print molbuild_residue_numbers
        dihedral_data = np.zeros(len(final_solution), dtype='int32, float32, float32')
        # print dihedral_data
        for i in xrange (0, len(dihedral_data)):
            dihedral_data[i][0] = int(molbuild_residue_numbers[i])
            dihedral_data[i][1] = final_solution[i][0]
            dihedral_data[i][2] = final_solution[i][1]
        
        rc.buildMolBuildInput(input_structure_path=molbuild_input, output_structure_path=molbuild_output, output_path='molbuild.in', mode=rc.MOL_BUILD_MODE_DIHEDRAL_ONLY, dihedral_data=dihedral_data)
        
           
def optimisePsiPhi(probabilityDistributions, probabilityDistributionPointers, solution, numIterations, move_lows, move_highs):
    solution_copy = np.copy(solution)
    # probs are in fortan order, psi, phi and not phi, psi
    for i in xrange (0, numIterations):
        for j in xrange (0, len(solution_copy)):
            
            currPhi = solution_copy[j][0]
            currPsi = solution_copy[j][1]
            
            
            probabilityDistribution = probabilityDistributions[probabilityDistributionPointers[j]]
            
            currProb = probabilityDistribution[int(np.floor((180 + currPsi) / 2))][int(np.floor((180 + currPhi) / 2))]
            psi_change = np.random.uniform(low=move_lows[j], high=move_highs[j])
            phi_change = np.random.uniform(low=move_lows[j], high=move_highs[j])

            new_phi = currPhi + phi_change
            new_psi = currPsi + psi_change

            if (new_psi < -180):
                new_psi = new_psi + 360

            if (new_psi >= 180):
                new_psi = new_psi - 360
                
            if (new_phi < -180):
                new_phi = new_phi + 360
                
            if (new_phi >= 180):
                new_phi = new_phi - 360
                
            grid_phi = int(np.floor((180 + new_phi) / 2))
            grid_psi = int(np.floor((180 + new_psi) / 2))
            if (currProb < probabilityDistribution[grid_psi][grid_phi] or np.random.rand() < probabilityDistribution[grid_psi][grid_phi]):
                solution_copy[j][0] = new_phi
                solution_copy[j][1] = new_psi
    return solution_copy


if __name__ == '__main__':
    doMain()
