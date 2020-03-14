'''
Replacers are used to modify a population after crossover and mutation have occured

.. codeauthor:: Nicholas Browning
'''

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import numpy as np
from . import Individual
import pygmo as pg
from src import outprocesses as op

def replace(settings, parents, children, replace_op, **args):
    '''

    Replace operation that selects the coressponding replacer function from the setting object

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    parents: list
        list of  :class:`src.gaapi.Individual.Individual`
    children: list
        list of  :class:`src.gaapi.Individual.Individual`
    replace_op: function
        the chosen replacer function in the :func:`main.initalise_ga`
    args

    Returns
    -------
    replace_op : function
        return values of the selected replacer

    '''
    return replace_op(settings, parents, children, **args)
    
def swap_replace(settings, parents, children, **args):
    '''

    Only returns the children

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    parents: list
        list of see :class:`src.gaapi.Individual.Individual`
    children: list
        list of see :class:`src.gaapi.Individual.Individual`
    args

    Returns
    -------

    '''
    return children

def elitist_replace(settings, parents, children, **args):
    '''

    Does an elitist replacement to improve GA convergence.
    According to the `elitist_factor` set a percentage of the most fit parents is retained.

    This is used only in single objective runs. The nsga2 algorithm is already elitist

    Use in the input file like so:

    .. code-block:: python

        [GA_REPLACER]
        replacer = elitism
        elitist_factor = 0.25

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    parents: list
        list of see :class:`src.gaapi.Individual.Individual`
    children: list
        list of see :class:`src.gaapi.Individual.Individual`
    args

    Returns
    -------
    elitst_population: list
        list of  :class:`src.gaapi.Individual.Individual` containing a population with elist parents chosen according to elistim_factor and the rest being children
    '''

    elitist_factor = settings.elitist_factor
    num_elites = int(elitist_factor * settings.population_size)

    print("ELITIST REPLACER with {} elite parents".format(num_elites))
    
    #parent_indexes = np.arange(settings.population_size)
    #child_indexes = np.arange(settings.population_size)  #deprecated
    
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


def non_dominated_sorting(settings, parents, children, **args):
    '''

    fast non dominated sorting (NDSA-II) by Deb et. al. to sort individuals.
    Automatically used if more than one evaluator is used.

    .. code-block:: python

        [GA_REPLACER]
        replacer = nsga2

    Parameters
    ----------
    settings : object
        see :class:`src.JobConfig.Settings`
    parents: list
        list of see :class:`src.gaapi.Individual.Individual`
    children: list
        list of see :class:`src.gaapi.Individual.Individual`
    args

    Returns
    -------
    new_population: list
        list containing the individuals for the next generation chosen according to non-domination rank and crowding distance sorting
    '''
    print("Using fast non dominated sorting (NDSA-II) by Deb et. al. to sort individuals")
    # make fitness array containing whole family
    familyFitnesses = [[0, 0]] * settings.population_size * 2  # TODO only two fitnesses currently. Need to make this general!
    for i in range(0, settings.population_size):
        familyFitnesses[i] = [parents[i].fitnesses[0], parents[i].fitnesses[1]]
        childi = i + settings.population_size
        familyFitnesses[childi] = [children[i].fitnesses[0], children[i].fitnesses[1]]
        pass
    if settings.verbose:
        print("Family Fitnesses", familyFitnesses)
        pass
    # sort using non dominated sorting and then select based on non-domination ranks and crowding distance, this returns the N best individuals that make up the next parent population
    new_population_idx = pg.select_best_N_mo(points=familyFitnesses, N=settings.population_size)
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(points=familyFitnesses)
    new_population = []
    if settings.verbose:
        print('id <', settings.population_size, '= parent population,', 'id >', settings.population_size, '= child population')
        print('nsga-II fittest individuals', new_population_idx)
        print('first front from nsga-II:')
    first_front = []
    ff_str = ''
    for x in ndf[0]:
        first_front.append(familyFitnesses[x])
        ff_str += str(familyFitnesses[x])
        if settings.store_intermediate_mol:
            if x > settings.population_size - 1:
                child_idx = x - settings.population_size
                molfile = children[int(child_idx)].mol[0]
                pass
            else:
                molfile = parents[int(x)].mol[0]
                pass
            indFitness = ' '
            for k in familyFitnesses[x]:
                indFitness += str(k) + ' '
                pass
            op.writeFrontFile(path=settings.output_path, file=molfile, fitness=indFitness,
                              iteration=settings.curr_iteration, frontid=x)
            pass
        pass

        pass
    if settings.verbose:
        print(ff_str)
        pass

    op.writeFrontLog(path=settings.output_path, front=first_front, iteration=settings.curr_iteration)
    pass

    # Loop through all individuals and include them in the new_population
    for i in new_population_idx:
        # because of the elitism we sort both the parents and children, we then have to add from the correct subpopulation (children or parents)
        if i > settings.population_size - 1:
            child_idx = i - settings.population_size
            new_population.append(children[int(child_idx)])
            pass
        else:
            new_population.append(parents[int(i)])

    return new_population
        
    
