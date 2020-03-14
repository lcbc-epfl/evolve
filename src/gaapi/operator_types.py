'''
@author: Nicholas Browning
'''
# #SUPPORTED OPERATOR TYPES#
import evaluators as evals
# terrible version of Enum struct as python 2.7 doesn't have it

# crossover types
SBX = 'sbx'
UNI_X = "uniform_crossover"

# mutator types
POLY_M = 'polynomial_mutation'
UNI_M = 'uniform_mutation'

# replacer types
SWAP = 'swap'
ELITIST = 'elitism'
NSGA2 = 'nsga2'

# generators
UNIFORM= 'uniform'
UNBIASED= 'random'
SWISSPROT= 'swissprot'

# selector ty[es
TOURNAMENT_WOR = 'tournament_wor'

availableEvaluators = {"helical_stability": evals.helical_stability,
                       "mmpbsa_multi": evals.mmpbsa_multi,
                       "stability_multi": evals.stability_multi,
                       "pmemd": evals.amber_energy_simplified}

