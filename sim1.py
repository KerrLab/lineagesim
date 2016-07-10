#!/usr/bin/env python3

import csv
import itertools
import json

import igraph
import sys
import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import multinomial as nmultinom
from numpy.random import normal as nnormal
from numpy.random import exponential as nexp
from six.moves import range as srange

from matplotlib import pyplot 


import cProfile, pstats
from io import StringIO
import pickle



POPSIZE = int(1e7)
NUM_CYCLES = 300
MUTATION_RATE = 1e-7
OUTFILENAME = "results.csv"
THRESH_FREQ = 0.001

np.random.seed(90210)

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
                     vertex_attrs = {'genotype_node': None,
                                     'first_seen': None,
                                     'last_seen': None,
                                     'depth': None,
                                     'abundance': 0,
                                     'abundances': [0],
                                     'frequency': 0,
                                     'frequencies': [],
                                     'max_frequency': 0,
                                     'fitness': 0,
                                     'mutation_node' : None,
                                     'fitness_diff': [0],
                                     'w_over_wbars' : [1]},
                     edge_attrs = {'relational':True})

    p.add_vertex(name = next(counter), depth = 0, abundance = population_size,
                 abundances = [population_size], frequency = 1.0,
                 max_frequency = 1.0, fitness = base_fitness,
                 fitness_diff = [0], w_over_wbars = [1], genotype_node = True, first_seen=0)
    return p


def reproduce(p, population_size):
    """Fill the population using fitness-proportional selection"""

    genotype_nodes = p.vs.select(genotype_node_eq=True)

    #first update abundance of genotypes
    fitnesses = np.array(genotype_nodes['fitness'])
    abundances = np.array(genotype_nodes['abundance'])
    ab_fit = fitnesses * abundances

    genotype_nodes['abundance'] = nmultinom(n = population_size,
                                  pvals = ab_fit / ab_fit.sum(),
                                  size = 1)[0]


    return p


def mutate_bdc(p, mutation_rate, node_counter, gen):
    """Mutate individuals in the population"""

    assert(mutation_rate >= 0 and mutation_rate <= 1)
    genotype_nodes = p.vs.select(genotype_node_eq=True)

    actual_pop_size = sum(genotype_nodes['abundance'])

    num_mutants = nbinom(n = genotype_nodes['abundance'], p = mutation_rate)
    genotype_nodes['abundance'] = genotype_nodes['abundance'] - num_mutants

    for parent_id in np.nonzero(num_mutants)[0]:
        for _mutant in srange(num_mutants[parent_id]):

            # TODO: handle mutation effect sizes properly
            #mu_effect = nnormal(loc = 0.0, scale = 0.1)

            mu_effect = nexp(scale=0.1)
            #add mutation to p
            new_mut_id = next(node_counter)

            p.add_vertex(name=new_mut_id,
                         genotype_node = False,
                         fitness = mu_effect,
                         frequency=0,
                         abundances=[],
                         max_frequency=0,
                         first_seen=gen)

            new_genotype_id = next(node_counter)

            p.add_vertex(name = new_genotype_id,
                         abundance = 1, abundances = [1],
                         depth = genotype_nodes[parent_id]['depth'] + 1,
                         fitness = genotype_nodes[parent_id]['fitness'] + mu_effect,
                         fitness_diff = [mu_effect], frequency = 0,
                         max_frequency = 0, genotype_node = True, first_seen=gen)
            
            #add edge to new mutation
            p.add_edge(source = p.vcount() - 2, target = p.vcount() - 1, relational=False)

            #add lineage edge between parent and offspring
            p.add_edge(source = genotype_nodes[parent_id], target = p.vcount() - 1,
                       fitness_effect = mu_effect, relational=True)

            #add mutational neighbors from parent to offspring vertex
            mutational_neighbors = p.vs[p.neighbors(genotype_nodes[parent_id])].select(genotype_node_eq=False)
            p.add_edges(zip(mutational_neighbors, [p.vcount()-1]*len(mutational_neighbors)))
    return p

def dilute(p, dilution_prob):
    """Thin the population"""
    genotype_nodes = p.vs.select(genotype_node_eq=True)

    assert(dilution_prob >= 0 and dilution_prob <= 1)
    genotype_nodes['abundance'] = nbinom(n = genotype_nodes['abundance'], p = dilution_prob)

    return p


def prune_frequency(p, min_frequency):
    """Delete extinct leaf nodes that did not reach a given frequency"""
    genotype_nodes = p.vs.select(genotype_node_eq=True)
    prunable_nodes = genotype_nodes.select(lambda v: v.outdegree() == 0 and v['abundance'] == 0 and v['first_seen'] is not None)

    prunable_nodes.delete()
    return p

def prune_mutations(p, min_frequency):
    mutation_nodes = p.vs.select(genotype_node_eq=False)
    mutation_nodes.select(lambda v: v['frequency']==0 and v['first_seen'] is not None and v['max_frequency'] < min_frequency).delete()
    return p


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
                             fieldnames = ['Generation', 'Genotype', 'Depth', 'Fitness', 'Abundance', 'Frequency'])
    outfile.writeheader()

    node_counter = itertools.count(0)

    genotypes = create_population(POPSIZE,  node_counter)

    for gen in srange(num_generations):
        genotype_nodes = genotypes.vs.select(genotype_node_eq=True)

        genotype_nodes.select(lambda v: v['abundance'] > 0)['last_seen'] = gen

        #print("Generation {c}. Max depth: {d}".format(c = gen, d = max(genotypes.vs['depth'])))
        #print("Gen {g}".format(g = gen))
        sys.stdout.write("\r")
        pct = float(gen) / num_generations
        sys.stdout.write("[%-40s] %d%%" % ('='* int(40 * pct), (pct * 100)))
        sys.stdout.flush()

        for extant in genotype_nodes.select(lambda v: v['abundance'] > 0):
            outfile.writerow({'Generation': gen,
                              'Genotype': extant.index,
                              'Depth': extant['depth'],
                              'Fitness': extant['fitness'],
                              'Abundance': extant['abundance'],
                              'Frequency': extant['frequency']})

        reproduce(genotypes, population_size = POPSIZE)
        mutate_bdc(genotypes, MUTATION_RATE, node_counter, gen)
        prune_frequency(genotypes, min_frequency = THRESH_FREQ)

        #update genotype subset after mutations...
        #TODO: is this required?
        genotype_nodes = genotypes.vs.select(genotype_node_eq=True)

        for v in genotype_nodes.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None):
            v['abundances'].append(int(v['abundance']))

        genotypes['generations'] += 1
        genotypes['population_size'] = int(sum(genotype_nodes['abundance']))
        genotype_nodes['frequency'] = genotype_nodes['abundance'] / sum(genotype_nodes['abundance'])

        for v in genotype_nodes.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']


        mutation_nodes = genotypes.vs.select(genotype_node_eq=False)
        #mutation_nodes.select(lambda v: v['first_seen'] is None)['first_seen'] = gen

        for m_node in mutation_nodes:
            m_node['frequency'] = np.sum(genotypes.vs[genotypes.neighbors(m_node)]['frequency'])
            #abusing the verticies def.
            m_node['abundances'].append(m_node['frequency'])

        for v in mutation_nodes.select(lambda v: v['frequency'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']


        prune_mutations(genotypes, min_frequency = THRESH_FREQ)
        
        #TODO: This line causes inconsistences! Bug! Ahhh! 
        prune_frequency(genotypes, min_frequency = THRESH_FREQ)

        #reset mutation nodes after pruning
        mutation_nodes = genotypes.vs.select(genotype_node_eq=False)

        mutation_nodes.select(lambda v: v['frequency'] > 0)['last_seen'] = gen


    #graph_write_json(genotypes, "tree-end.json", sort_keys = True)
    genotypes.write_gml("tree-end.gml")
    pickle.dump(genotypes, open('tree-end.pickle', 'wb'))
    return(genotypes)
# -----------------------------------------------------------------------------

'''
if __name__ == "__main__":
    run_simulation(num_generations = NUM_CYCLES)
'''

pr = cProfile.Profile()

pr.enable()
# ... do something ...

genotypes = run_simulation(num_generations = NUM_CYCLES)
mutation_nodes = genotypes.vs.select(genotype_node_eq = False)
print (mutation_nodes)


ms = (mutation_nodes['abundances'])
mt = (mutation_nodes['first_seen'])

pr.disable()

xvec = range(0, NUM_CYCLES)
yvec = [0 for i in range(0, NUM_CYCLES)]
for i_mutant in range(len(ms)):
    yvec = [0 for i in range(0, NUM_CYCLES)]
    abundances = ms[i_mutant]
    first_seen = mt[i_mutant]
    for j in range(len(abundances)):
        yvec[j+first_seen] = abundances[j]

    print(mutation_nodes[i_mutant])



