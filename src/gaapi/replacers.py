'''
crossovers.py

@author: Nicholas Browning
'''
import numpy as np
import Individual

def replace(settings, parents, children, replace_op, **args):
    return replace_op(settings, parents, children, **args)
    

def swap_replace(settings, parents, children, **args):
    return children

def elitist_replace(settings, parents, children, **args):
    
    elitist_factor = settings.elitist_factor
    
    num_elites = int(elitist_factor * settings.population_size)
    
    parent_indexes = np.arange(settings.population_size)
    child_indexes = np.arange(settings.population_size)
    
    sorted_parent_indexes = np.argsort(parents)

    sorted_child_indexes = np.argsort(children)
    
    if (settings.verbose):
        print "elite parent indexes:", sorted_parent_indexes[:num_elites]
        print "elite children indexes:", sorted_child_indexes[:settings.population_size - num_elites]
    
    
    elitist_population = []
    
    for i in xrange (num_elites):
        indiv = Individual.Individual(settings, parents[sorted_parent_indexes[i]])
        elitist_population.append(indiv)
    
    for i in xrange (settings.population_size - num_elites):
        indiv = Individual.Individual(settings, children[sorted_child_indexes[i]])
        elitist_population.append(indiv)
        
        
    return elitist_population
        
    
