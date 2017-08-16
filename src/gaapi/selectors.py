'''
crossovers.py

@author: Nicholas Browning
'''

import numpy as np

def selector(settings, individuals, selector_op, **args):
    return selector_op(settings, individuals, **args)


def tournamentSelectionWOR(settings, individuals, **kargs):
    tournament_size = settings.tournament_size
    
    indexes = np.arange(settings.population_size)
    
    result = []
    
    for _ in xrange(0, settings.population_size):
        np.random.shuffle(indexes)
        min_idx = np.argmin([individuals[i] for i in indexes[:tournament_size]])
        result.append(indexes[:tournament_size][min_idx])
        
    return np.asarray(result)
