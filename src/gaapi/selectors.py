'''
This contains the selectors used for  selecting the individuals allowed to mate

.. codeauthor:: Nicholas Browning
'''
from __future__ import print_function

from builtins import range
import numpy as np

def selector(settings, individuals, selector_op, **args):
    '''

    Prints the selector used to the output and run the selector

    Parameters
    ----------
    settings: object
        see :class:`src.JobConfig.Settings`
    individuals: list
        list of individuals of this generation
    selector_op : str
        defined by JobInput
    args

    Returns
    -------
    selector_op : func
        selected selctor function

    '''
    print("SELECTORS {}".format(settings.selector))
    return selector_op(settings, individuals, **args)


def tournamentSelectionWOR(settings, individuals, **kargs):
    '''

    Tournament selection where two individuals fight a tournament and the stronger one is allowed to mate

    .. code-block:: python

        [GA_SELECTION]
        selector = tournament_wor
        tournament_size = 2

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    individuals
    kargs

    Returns
    -------
    indeces: np.array
        indeces of selected individuals

    '''
    tournament_size = settings.tournament_size
    
    indexes = np.arange(settings.population_size)
    
    result = []
    
    for _ in range(0, settings.population_size):
        np.random.shuffle(indexes)
        min_idx = np.argmin([individuals[i] for i in indexes[:tournament_size]])
        result.append(indexes[:tournament_size][min_idx])
        
    return np.asarray(result)

def roulette_wheel(settings, individuals, **kargs):
    '''

    not implemented

    Parameters
    ----------
    settings
    individuals
    kargs

    Returns
    -------

    '''
    pass
