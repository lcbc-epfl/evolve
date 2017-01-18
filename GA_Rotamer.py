'''
GA_Rotamer.py

@author: Nicholas Browning
'''

#!/software/anaconda-2.0.1/bin/python

'''
#!/software/anaconda-2.0.1/bin/python
import inspyred as insp
insp.ec.analysis.generation_plot(open("stats.csv", 'r'))
'''


from random import Random
from time import time
import inspyred as insp
import RotamerCreator as rc
import matplotlib.pyplot as plt
import numpy as np


'''
Created on Feb 2, 2016

V1.1

@author: Nicholas Browning
'''

fitnessListFile = "fitness.dat"
# CHROMOSOME BOUNDS DEFINED HERE

# GLY GLY PRO GLY
psiphi_nopro_data = [[26, 14, 15, 16], [26, 14, 13, 28], [24, 10, 11, 12], [24, 10, 9, 31], [34, 33, 1, 32], [34, 33, 36, 39]]
psiphi_data = [[26, 14, 15, 16], [26, 14, 13, 28], [24, 10, 11, 12], [24, 10, 9, 31], [6, 3, 4, 5], [6, 3, 2, 1], [34, 33, 1, 32], [34, 33, 36, 39]]

lowerBounds = [-180, -180, -180, -180, -180, -180, -180, -180]
upperBounds = [180, 180, 180, 180, 180, 180, 180, 180]

move_lows = [-25, -25, -25, -25, -100, -100, -50, -50]
move_highs = [25, 25, 25, 25, 100, 100, 50, 50]

prob_pointers = [1, 1, 1, 1, 3, 3, 2, 2]
POP_SIZE = 10
initial_structure = "TOP5_CIS.pdb"
fitnessListC = []
fitnessListE = []

prng = Random()
prng.seed(time()) 

ea = insp.ec.EvolutionaryComputation(prng)
 
def main(display=False): 

    rc.initRotamerLists(dielectric=10)

    ea.selector = insp.ec.selectors.tournament_selection
    ea.variator = [insp.ec.variators.simulated_binary_crossover, insp.ec.variators.polynomial_mutation]
    ea.replacer = insp.ec.replacers.generational_replacement
    ea.terminator = [insp.ec.terminators.generation_termination]
    
    ea.observer = [insp.ec.observers.stats_observer, insp.ec.observers.file_observer]
    # stats_observer prints generation statistics to terminal
    # file_observer prints stats and individuals to files: statistics_file=open("stats.csv", 'w'), individuals_file=open("indivs.csv", 'w')
   
    final_pop = ea.evolve(bounder=insp.ec.Bounder(lower_bound=lowerBounds, upper_bound=upperBounds), generator=generate, evaluator=evaluateFitness, maximize=False,
    max_generations=5, num_elites=5, pop_size=10,
    num_selected=10, tournament_size=2,
    crossover_rate=1.0, sbx_distribution_index=5, discrete_crossover=True,
    eta=5, mutation_rate=0.15, discrete_mutation=True, statistics_file=open("stats.csv", 'w'), individuals_file=open("indivs.csv", 'w')
    
    
    )
                          
    if display:
        best = max(final_pop)
        print('Best Solution: \n{0}'.format(str(best)))
    return ea


            
def defaultGenerate(random, args):
    return [random.randint(lowerBounds[i], upperBounds[i]) for i in xrange (0, len(upperBounds))]

def monteCarloGenerate(random, solution):
    import montecarlo.optimise as opt
    probs = opt.loadDefaultProbabilityDistributions()
    return opt.optimisePsiPhi(probs, prob_pointers, solution, 10000, move_lows, move_highs)

def generate(random, args):
    return monteCarloGenerate(random, defaultGenerate(random, args))

def evaluateFitness(candidates, args):
    import shutil as sh
    import os
    parent_indexes = args.setdefault('parent_indexes', None)
    fitnesses = []
    
    if (not os.path.exists("temp")):
        os.mkdir("temp")
        
    for i, c in enumerate (candidates):
        print "Evaluating candidate:", c
        
        init_structure = initial_structure
        
        if (parent_indexes is not None):
            if (parent_indexes[i]) == -1:
                # use initial structure
                pass
            elif (parent_indexes[i] > POP_SIZE):
                init_structure = "opt_child_" + str(ea.num_generations - 2) + "_" + str(parent_indexes[i] - POP_SIZE) + ".xyz"
            else:
                init_structure = "opt_child_" + str(ea.num_generations - 1) + "_" + str(parent_indexes[i]) + ".xyz"
        index = -1
        if (fitnessListC is not None):
            try:
                index = fitnessListC.index(c)
            except ValueError:
                pass
        if (index == -1):
            name = "child" + "_" + str(ea.num_generations) + "_" + str(i)
            # fitness = rc.evaluateRotamerLEAPGeneration(initialMoleculePath=initialMolPath, originalEnergy=origE, moleculeSubstitutionIndexes=molindexes, designVector=c)
            # fitness = rc.evaluateRotamerMolBuilderGeneration(initialMoleculePath=initialMolPath, originalEnergy=origE, molecule_start_end_indexes=mol_alphabetacarbonlist, designVector=c)

            fitness = rc.evaluateRotamerMolBuilderGaussian(initialMoleculePath=init_structure, gaussianOutputName=name, outputMolbuildFile=name + ".molbuild", outputStructurePath=name + '.xyz' , designVector=c, m=rc.MOL_BUILD_MODE_DIHEDRAL_ONLY, data=psiphi_data)
            fitnesses.append(fitness)
            rc.extractOptimisedGeometry(logfile=name + '.log', out_structure='opt_' + name + '.xyz')
            # print getBackboneDihedrals(dihedral_data=psiphi_data, structure_in='opt_solution.xyz', molbuild_in='test.mb') 
            if (fitnessListC is not None):
                fitnessListC.append(c)
                fitnessListE.append(fitness)
                fitnessListOutput = open(fitnessListFile, 'a')
                fitnessListOutput.write(" ".join(str(d) for d in c) + ", " + str(fitness) + "\n")
                fitnessListOutput.close()
                
        else:
            fitnesses.append(fitnessListE[index])
            print "Updating solution:", c , " with fitness:", fitnessListE[index]
            
    return fitnesses

def testDisplay():
    from montecarlo import optimise
    probs = optimise.loadDefaultProbabilityDistributions()
 
    sols = np.zeros((50, 8))
    mc_sols = np.zeros((50, 8))    
     
    for i in xrange (0, 50):
        sols[i] = defaultGenerate(prng, None)
        mc_sols[i] = sols[i]
         
    for i in xrange (0, 50):
        mc_sols[i] = monteCarloGenerate(prng, sols[i])
   
    optimise.plotProbabilityDistributionAndSolution(probs, [1, 1, 1, 1, 3, 3, 2, 2], mc_sols)
       
    moleculedata, fragdata, dihedraldata = rc.createdata(sols[49], psiphi_data, rc.MOL_BUILD_MODE_DIHEDRAL_ONLY)
    rc.buildMolBuildInput(input_structure_path="TOP5_CIS.pdb", output_path="test.molbuild", mode=rc.MOL_BUILD_MODE_DIHEDRAL_ONLY, molecule_data=moleculedata, frag_data=fragdata, subst_options=None, dihedral_data=dihedraldata, output_structure_path='before.xyz')
    rc.runMolBuild(inputFile='test.molbuild')
     
    moleculedata, fragdata, dihedraldata = rc.createdata(mc_sols[49], psiphi_data, rc.MOL_BUILD_MODE_DIHEDRAL_ONLY)
    rc.buildMolBuildInput(input_structure_path="TOP5_CIS.pdb", output_path="test.molbuild", mode=rc.MOL_BUILD_MODE_DIHEDRAL_ONLY, molecule_data=moleculedata, frag_data=fragdata, subst_options=None, dihedral_data=dihedraldata, output_structure_path='after.xyz')
    rc.runMolBuild(inputFile='test.molbuild')
    
if __name__ == '__main__':
    print "Start"
    prng = Random()
    prng.seed(time()) 
    
    # uncomment the following with the insp.ec.observers.plot_observer to get the interactive plot
    # plt.ioff()
    # plt.show()
    
    # uncomment the following to get the allele plot and full convergence graph after simulation
    # insp.ec.analysis.allele_plot("indivs.csv")
    # insp.ec.analysis.generation_plot(open("stats.csv", 'r'))
    # testDisplay()
    main(display=True)
    print "Finish"
