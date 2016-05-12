#!/usr/bin/env python3

import csv
import itertools
import json
import sys

import igraph
import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import multinomial as nmultinom
from numpy.random import normal as nnormal
from six.moves import range as srange

from mutate_multiples import mutate_multiples

POPSIZE = int(1e6)
NUM_CYCLES = 1000
MUTATION_RATE = 1e-6
OUTFILENAME = "results.csv"
THRESH_FREQ = 0.001

#np.random.seed(90210)

# current abundances: [v['abundances'][-1] for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None)]
# g.vs(abundance_gt=0)
# g.vs[1].index

# TODO: use abundances[-1] instead of abundance to get current abundance?

# -----------------------------------------------------------------------------

# Create the population
def create_population(population_size, counter, base_fitness = 1.0):
    """Create the population"""
    p = igraph.Graph(directed = True,
                     graph_attrs = {'population_size': population_size,
                                    'generations': 0},
                     vertex_attrs = {'first_seen': None,
                                     'last_seen': None,
                                     'depth': None,
                                     'abundance': 0,
                                     'abundances': [0],
                                     'total_abundance': 0,
                                     'frequency': 0,
                                     'max_frequency': 0,
                                     'fitness': 0,
                                     'fitness_diff': [0]})

    p.add_vertex(name = next(counter), depth = 0, abundance = population_size,
                 abundances = [population_size],
                 total_abundance = population_size, frequency = 1.0,
                 max_frequency = 1.0, fitness = base_fitness,
                 fitness_diff = [0])
    return p


def reproduce(p, population_size):
    """Fill the population using fitness-proportional selection"""

    fitnesses = np.array(p.vs['fitness'])
    abundances = np.array(p.vs['abundance'])
    ab_fit = fitnesses * abundances

    p.vs['abundance'] = nmultinom(n = population_size,
                                  pvals = ab_fit / ab_fit.sum(),
                                  size = 1)[0]
    return p


def mutate_bdc(p, mutation_rate, genotype_counter):
    """Mutate individuals in the population"""

    assert(mutation_rate >= 0 and mutation_rate <= 1)

    num_mutants = nbinom(n = p.vs['abundance'], p = mutation_rate)
    p.vs['abundance'] = p.vs['abundance'] - num_mutants

    for parent_id in np.nonzero(num_mutants)[0]:
        for _mutant in srange(num_mutants[parent_id]):

            # TODO: handle mutation effect sizes properly
            mu_effect = nnormal(loc = 0.0, scale = 0.1)

            p.add_vertex(name = next(genotype_counter),
                         abundance = 1, abundances = [1],
                         total_abundance = 1,
                         depth = p.vs[parent_id]['depth'] + 1,
                         fitness = p.vs[parent_id]['fitness'] + mu_effect,
                         fitness_diff = [mu_effect], frequency = 0,
                         max_frequency = 0)
            p.add_edge(source = parent_id, target = p.vcount() - 1,
                       fitness_effect = mu_effect)

    return p


def dilute(p, dilution_prob):
    """Thin the population"""
    assert(dilution_prob >= 0 and dilution_prob <= 1)
    p.vs['abundance'] = nbinom(n = p.vs['abundance'], p = dilution_prob)
    return p


def prune_frequency(p, min_frequency):
    """Delete extinct leaf nodes that did not reach a given frequency"""
    p.vs.select(lambda v: v.outdegree() == 0 and v['abundance'] == 0 and v['first_seen'] is not None and v['max_frequency'] < min_frequency).delete()
    return p


def get_total_abundances(v):
    """Get the abundance of a vertex and all of its descendants"""

    v['total_abundance'] = v['abundance']

    for child in v.neighbors(mode = "OUT"):
        v['total_abundance'] += get_total_abundances(v = child)

    return v['total_abundance']

# -----------------------------------------------------------------------------

def graph_write_json(g, filename, **kwargs):
    """Write the graph to JSON file"""
    d = {}
    d['attributes'] = {a:g[a] for a in g.attributes()}
#    d['vertices'] = [{'index': v.index,
#                      'attributes': v.attributes()} for v in g.vs]
    d['vertices'] = [{'index': v.index,
                      'attributes': {'name': v['name'],
                                     'first_seen': v['first_seen'],
                                     'last_seen': v['last_seen'],
                                     'depth': v['depth'],
                                     'abundance': int(v['abundance']),
                                     'abundances': v['abundances'],
                                     'total_abundance': int(v['total_abundance']),
                                     'frequency': v['frequency'],
                                     'max_frequency': v['max_frequency'],
                                     'fitness': v['fitness'],
                                     'fitness_diff': v['fitness_diff']}} for v in g.vs]

    d['edges'] = [{'index': e.index,
                   'source': e.source,
                   'target': e.target,
                   'attributes': e.attributes()} for e in g.es]

    with open(filename, 'w') as outfile:
        json.dump(d, outfile, **kwargs)

# -----------------------------------------------------------------------------

def run_simulation(num_generations):

    outfile = csv.DictWriter(open(OUTFILENAME, 'w'),
                             fieldnames = ['Generation', 'Genotype', 'Depth', 'Fitness', 'Abundance', 'TotAbundance', 'Frequency'])
    outfile.writeheader()

    genotype_counter = itertools.count(0)
    genotypes = create_population(population_size = POPSIZE, counter = genotype_counter)

    for gen in srange(num_generations):
        genotypes.vs.select(lambda v: v['first_seen'] is None)['first_seen'] = gen
        genotypes.vs.select(lambda v: v['abundance'] > 0)['last_seen'] = gen

        #print("Generation {c}. Max depth: {d}".format(c = gen, d = max(genotypes.vs['depth'])))
        #print("Gen {g}".format(g = gen))
        sys.stdout.write("\r")
        pct = float(gen) / num_generations
        sys.stdout.write("[%-40s] %d%%" % ('='* int(40 * pct), (pct * 100)))
        sys.stdout.flush()

        for g in genotypes.vs.select(lambda v: v['total_abundance'] > 0):
            outfile.writerow({'Generation': gen,
                              'Genotype': g.index,
                              'Depth': g['depth'],
                              'Fitness': g['fitness'],
                              'Abundance': g['abundance'],
                              'TotAbundance': g['total_abundance'],
                              'Frequency': g['frequency']})

        reproduce(genotypes, population_size = POPSIZE)
        #mutate_multiples(genotypes, mutation_rate = MUTATION_RATE, genotype_counter = genotype_counter)
        mutate_bdc(genotypes, mutation_rate = MUTATION_RATE, genotype_counter = genotype_counter)
        prune_frequency(genotypes, min_frequency = THRESH_FREQ)
        get_total_abundances(v = genotypes.vs[0])

        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None):
            v['abundances'].append(int(v['abundance']))

        genotypes['generations'] += 1
        genotypes['population_size'] = int(sum(genotypes.vs['abundance']))
        genotypes.vs['frequency'] = genotypes.vs['abundance'] / sum(genotypes.vs['abundance'])

        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']

        # For writing the tree at every cycke
        #genotypes.write_gml("TREES/genotypes-{0:06d}.gml".format(gen))

    graph_write_json(genotypes, "tree-end.json", sort_keys = True)
    genotypes.write_gml("tree-end.gml")

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    run_simulation(num_generations = NUM_CYCLES)

