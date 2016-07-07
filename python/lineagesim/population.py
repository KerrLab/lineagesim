# -*- coding: utf-8 -*-

"""Functions for creating and manipulating populations of genotypes"""

import igraph

import numpy as np
from numpy import hstack as nhstack
from numpy import nonzero as nnonzero
from numpy import prod as nprod
from numpy import subtract as nsubtract
from numpy.random import binomial as nbinom
from numpy.random import exponential as nexp
from numpy.random import multinomial as nmultinom

from scipy.special import gammainc as gamma

from six.moves import range as srange


def create_population(population_size, counter, base_fitness=1.0):
    """Create the population"""
    pop = igraph.Graph(directed=True,
                       graph_attrs={'population_size': population_size,
                                    'generations': 0},
                       vertex_attrs={'name': None,
                                     'parent': None,
                                     'first_seen': None,
                                     'last_seen': None,
                                     'lineage_last_seen': None,
                                     'depth': None,
                                     'abundance': 0,
                                     'abundances': [0],
                                     'lineage_abundance': 0,
                                     'frequency': 0,
                                     'max_frequency': 0,
                                     'lineage_frequency': 0,
                                     'fixation_time': -1,
                                     'fitness': 0,
                                     'fitness_effects': [0],
                                     'fitness_diff': 0})

    pop.add_vertex(name=next(counter),
                   parent=-1,
                   depth=0,
                   abundance=population_size,
                   abundances=[population_size],
                   lineage_abundance=population_size,
                   frequency=1.0,
                   max_frequency=1.0,
                   lineage_frequency=1.0,
                   fixation_time=0,
                   fitness=base_fitness,
                   fitness_effects=[0],
                   fitness_diff=0)
    return pop


def reproduce(population, population_size):
    """Fill the population using fitness-proportional selection"""

    fitnesses = np.array(population.vs['fitness'])
    abundances = np.array(population.vs['abundance'])
    ab_fit = fitnesses * abundances

    population.vs['abundance'] = nmultinom(n=population_size,
                                           pvals=ab_fit / ab_fit.sum(),
                                           size=1)[0]
    return population



def mutate(population, mutation_rate, mutation_scale, genotype_counter):
    """Mutate individuals in the population"""

    assert mutation_rate >= 0 and mutation_rate <= 1

    num_mutants = nbinom(n=population.vs['abundance'],
                         p=1 - np.exp(-mutation_rate))
    population.vs['abundance'] = population.vs['abundance'] - num_mutants

    r = num_mutants
    k = 1
    while any(r > 0):
        # r is the number of individuals with at least k+1 mutations
        # rate inside binomial is the prob of an individual having at least
        # k+1 mutations given that it has at least k mutations
        new_r = nbinom(n=r,
                       p=gamma(k+1, mutation_rate)/gamma(k, mutation_rate),
                       size=len(r))
        num_k_mutants = nsubtract(r, new_r) #number of k-mutants
        for parent_id in nnonzero(num_k_mutants)[0]:
            for _mutant in srange(num_k_mutants[parent_id]):
                mu_effect = nexp(scale=mutation_scale, size=k)
                mutant_fitness = nhstack((population.vs[parent_id]['fitness'], (1.0 + mu_effect))).prod()

                population.add_vertex(name=next(genotype_counter),
                                      parent=int(population.vs[parent_id]['name']),
                                      abundance=1,
                                      abundances=[1],
                                      total_abundance=1,
                                      depth=population.vs[parent_id]['depth'] + 1,
                                      fitness=mutant_fitness,
                                      fitness_effects=mu_effect.tolist(),
                                      fitness_diff=mutant_fitness - population.vs[parent_id]['fitness'],
                                      frequency=1.0 / population['population_size'],
                                      max_frequency=1.0 / population['population_size'],
                                      fixation_time=-1)
                population.add_edge(source=parent_id,
                                    target=population.vcount() - 1,
                                    fitness_effect=nprod(1 + mu_effect))

        r = new_r
        k += 1
    return population


def prune_frequency(population, min_frequency):
    """Delete extinct leaf nodes that did not reach a given frequency"""
    # Could also prune based on lineage_abundance and max_frequency
    population.vs.select(lambda v: v.outdegree() == 0 and v['abundance'] == 0 and v['first_seen'] is not None and v['max_frequency'] < min_frequency).delete()
    return population


def get_lineage_abundances(genotype):
    """Get the abundance of a genotype and all of its descendants"""

    genotype['lineage_abundance'] = genotype['abundance']

    for child in genotype.neighbors(mode="OUT"):
        genotype['lineage_abundance'] += get_lineage_abundances(genotype=child)

    return genotype['lineage_abundance']

