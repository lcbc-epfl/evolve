'''
Replacers are used to modify a population after crossover and mutation have occured

.. codeauthor:: Nicholas Browning
'''

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import numpy as np
from . import Individual

def replace(settings, parents, children, replace_op, **args):
    '''

    Replace operation that selects the coressponding replacer function from the setting object

    Parameters
    ----------
    settings
    parents
    children
    replace_op
    args

    Returns
    -------
    replace_op : function
        selected replacer

    '''
    return replace_op(settings, parents, children, **args)
    
def swap_replace(settings, parents, children, **args):
    '''

    Only returns the children

    Parameters
    ----------
    settings
    parents
    children
    args

    Returns
    -------

    '''
    return children

def elitist_replace(settings, parents, children, **args):
    '''

    Does an elitist replacement to improve GA convergence.
    According to the `elitist_factor` set a percentage of the most fit parents is retained.

    Use in the input file like so:

        [GA_REPLACER]
        replacer = elitism
        elitist_factor = 0.25

    Parameters
    ----------
    settings
    parents
    children
    args

    Returns
    -------

    '''

    elitist_factor = settings.elitist_factor
    num_elites = int(elitist_factor * settings.population_size)

    print("ELITIST REPLACER with {} elite parents".format(num_elites))
    
    parent_indexes = np.arange(settings.population_size)
    child_indexes = np.arange(settings.population_size)
    
    sorted_parent_indexes = np.argsort(parents)
    sorted_child_indexes = np.argsort(children)
    
    if (settings.verbose):
        print("elite parent indexes:", sorted_parent_indexes[:num_elites])
        print("elite children indexes:", sorted_child_indexes[:settings.population_size - num_elites])
    
    
    elitist_population = []
    
    for i in range (num_elites):
        indiv = Individual.Individual(settings, parents[sorted_parent_indexes[i]])
        elitist_population.append(indiv)
    
    for i in range (settings.population_size - num_elites):
        indiv = Individual.Individual(settings, children[sorted_child_indexes[i]])
        elitist_population.append(indiv)
        
        
    return elitist_population
        
    
