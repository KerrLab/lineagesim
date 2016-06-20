#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "0.1.0"

import argparse
import datetime
import csv
import itertools
import json
import os
import sys
from warnings import warn

import igraph
import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import multinomial as nmultinom
from numpy.random import normal as nnormal
from six.moves import range as srange

from mutate_multiples import mutate_multiples

OUTFILENAME = "results.csv"

# current abundances: [v['abundances'][-1] for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None)]
# g.vs(abundance_gt=0)
# g.vs[1].index

# TODO: use abundances[-1] instead of abundance to get current abundance?

# -----------------------------------------------------------------------------

# Create the population
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
                                     'total_abundance': 0,
                                     'frequency': 0,
                                     'max_frequency': 0,
                                     'total_frequency': 0,
                                     'fixation_time': -1,
                                     'fitness': 0,
                                     'fitness_effects': [0],
                                     'fitness_diff': 0})

    pop.add_vertex(name=next(counter),
                   parent=-1,
                   depth=0,
                   abundance=population_size,
                   abundances=[population_size],
                   total_abundance=population_size,
                   frequency=1.0,
                   max_frequency=1.0,
                   total_frequency=1.0,
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


def mutate_bdc(population, mutation_rate, genotype_counter):
    """Mutate individuals in the population"""

    assert mutation_rate >= 0 and mutation_rate <= 1

    num_mutants = nbinom(n=population.vs['abundance'], p=mutation_rate)
    population.vs['abundance'] = population.vs['abundance'] - num_mutants

    for parent_id in np.nonzero(num_mutants)[0]:
        for _mutant in srange(num_mutants[parent_id]):

            # TODO: handle mutation effect sizes properly
            mu_effect = nnormal(loc=0.0, scale=0.1)

            population.add_vertex(name=next(genotype_counter),
                                  parent=int(population.vs[parent_id]['name']),
                                  abundance=1,
                                  abundances=[1],
                                  total_abundance=1,
                                  depth=population.vs[parent_id]['depth'] + 1,
                                  fitness=population.vs[parent_id]['fitness'] + mu_effect,
                                  fitness_effects=[mu_effect],
                                  fitness_diff=mu_effect,
                                  frequency=1 / population['population_size'],
                                  fixation_time=-1,
                                  max_frequency=0)
            population.add_edge(source=parent_id,
                                target=population.vcount() - 1,
                                fitness_effect=mu_effect)

    return population


def dilute(population, dilution_prob):
    """Thin the population"""
    assert dilution_prob >= 0 and dilution_prob <= 1
    population.vs['abundance'] = nbinom(n=population.vs['abundance'],
                                        p=dilution_prob)
    return population


def prune_frequency(population, min_frequency):
    """Delete extinct leaf nodes that did not reach a given frequency"""
    # Could also prune based on total_abundance and max_frequency
    population.vs.select(lambda v: v.outdegree() == 0 and v['abundance'] == 0 and v['first_seen'] is not None and v['max_frequency'] < min_frequency).delete()
    return population


def get_total_abundances(genotype):
    """Get the abundance of a genotype and all of its descendants"""

    genotype['total_abundance'] = genotype['abundance']

    for child in genotype.neighbors(mode="OUT"):
        genotype['total_abundance'] += get_total_abundances(genotype=child)

    return genotype['total_abundance']

# -----------------------------------------------------------------------------

def graph_write_json(graph, filename, **kwargs):
    """Write the graph to JSON file"""
    d = {}
    d['attributes'] = {a:graph[a] for a in graph.attributes()}
#    d['vertices'] = [{'index': v.index,
#                      'attributes': v.attributes()} for v in graph.vs]
    d['vertices'] = [{'index': v.index,
                      'attributes': {'name': v['name'],
                                     'parent': v['parent'],
                                     'first_seen': v['first_seen'],
                                     'last_seen': v['last_seen'],
                                     'lineage_last_seen': v['lineage_last_seen'],
                                     'depth': v['depth'],
                                     'abundance': int(v['abundance']),
                                     'abundances': v['abundances'],
                                     'total_abundance': int(v['total_abundance']),
                                     'frequency': v['frequency'],
                                     'max_frequency': v['max_frequency'],
                                     'total_frequency': v['total_frequency'],
                                     'fixation_time': v['fixation_time'],
                                     'fitness': v['fitness'],
                                     'fitness_effects': v['fitness_effects'],
                                     'fitness_diff': v['fitness_diff']}} for v in graph.vs]

    d['edges'] = [{'index': e.index,
                   'source': e.source,
                   'target': e.target,
                   'attributes': e.attributes()} for e in graph.es]

    with open(filename, 'w') as outfile:
        json.dump(d, outfile, **kwargs)


def graph_write_csv(graph, filename, **kwargs):
    """Write the graph to CSV file"""

    with open(filename, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['Genotype',
                                                     'Parent',
                                                     'FirstSeen',
                                                     'LastSeen',
                                                     'LineageLastSeen',
                                                     'Depth',
                                                     'Fitness',
                                                     'FitnessEffects',
                                                     'FitnessDiff',
                                                     'Abundance',
                                                     'TotalAbundance',
                                                     'Frequency',
                                                     'MaxFrequency',
                                                     'TotalFrequency',
                                                     'FixationTime'])
        writer.writeheader()

        for g in graph.vs.select(lambda v: v['first_seen'] is not None):
            writer.writerow({'Genotype': g['name'],
                             'Parent': g['parent'],
                             'FirstSeen': g['first_seen'],
                             'LastSeen': g['last_seen'],
                             'LineageLastSeen': g['lineage_last_seen'],
                             'Depth': g['depth'],
                             'Fitness': g['fitness'],
                             'FitnessEffects': sum(g['fitness_effects']),
                             'FitnessDiff': g['fitness_diff'],
                             'Abundance': g['abundance'],
                             'TotalAbundance': g['total_abundance'],
                             'Frequency': g['frequency'],
                             'MaxFrequency': g['max_frequency'],
                             'TotalFrequency': g['total_frequency'],
                             'FixationTime': g['fixation_time']})


# -----------------------------------------------------------------------------

def parse_arguments():
    """Parse command line arguments"""

    def check_positive_01(value):
        """Make sure a value is between 0 and 1, inclusive"""
        fval = float(value)
        if fval < 0 or fval > 1:
            raise argparse.ArgumentTypeError("%s is not in the range [0,1]" % value)
        return fval

    def check_nonnegative(value):
        """Make sure a value is non-negative"""
        fval = float(value)
        if fval < 0:
            raise argparse.ArgumentTypeError("%s is not a valid non-negative value" % value)
        return fval

    def check_nonnegative_int(value):
        """Make sure an integer value is non-negative"""
        ival = int(value)
        if ival < 0:
            raise argparse.ArgumentTypeError("%s is not a valid non-negative integer value" % value)
        return ival

    def check_positive_int(value):
        """Make sure an integer value is positive"""
        ival = int(value)
        if ival <= 0:
            raise argparse.ArgumentTypeError("%s is not a valid positive integer value" % value)
        return ival


    parser = argparse.ArgumentParser(prog='lineagesim',
                                     description='Run a lineage interference simulation')
    parser.add_argument('--data_dir', '-d', metavar='DIR', default='data',
                        help='Directory to store data (default: data)')
    parser.add_argument('--fixation_freq', '-f', metavar='F',
                        type=check_positive_01,
                        help='Threshold frequency for classification as fixed (default: 1 - μ)')
    parser.add_argument('--generations', '-G', metavar='G', default=1000,
                        type=check_positive_int,
                        help='Number of generations to simulate')
    parser.add_argument('--mutation_rate', '-m', metavar='μ', default=1e-6,
                        type=check_positive_01,
                        help='Mutation rate (default: 1e-6)')
    parser.add_argument('--population_size', '-N', metavar='N',
                        default=int(1e6), type=check_positive_int,
                        help='Size of the population (default: 1,000,000)')
    parser.add_argument('--prune_freq', '-p', metavar='F',
                        type=check_positive_01,
                        help='Threshold frequency for pruning (default: (N * μ) / N)')
    parser.add_argument('--seed', '-s', metavar='S', help='Set the '\
                        'pseudorandom number generator seed',
                        type=check_nonnegative_int)
    parser.add_argument('--quiet', '-q', action='store_true', default=False,
                        help='Suppress output messages')
    parser.add_argument('--version', action='version', version=__version__)

    return parser.parse_args()


def run_simulation(args=parse_arguments()):
    """Run the simulation"""

    if not args.fixation_freq:
        args.fixation_freq = 1.0 - args.mutation_rate
    if not args.prune_freq:
        args.prune_freq = (args.population_size * args.mutation_rate) / args.population_size

    if args.seed:
        np.random.seed(args.seed)

    if os.path.exists(args.data_dir):
        newname = '{o}-{d}'.format(o=args.data_dir, d=datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
        msg = "{d} already exists. Saving data to {new}.".format(d=args.data_dir, new=newname)
        warn(msg)
        args.data_dir = newname

    os.mkdir(args.data_dir)


    outfile = csv.DictWriter(open(os.path.join(args.data_dir, OUTFILENAME),
                                  'w'),
                             fieldnames=['Generation',
                                         'Genotype',
                                         'Parent',
                                         'FirstSeen',
                                         'LastSeen',
                                         'LineageLastSeen',
                                         'Depth',
                                         'Fitness',
                                         'FitnessEffects',
                                         'FitnessDiff',
                                         'Abundance',
                                         'TotalAbundance',
                                         'Frequency',
                                         'MaxFrequency',
                                         'TotalFrequency',
                                         'FixationTime'])
    outfile.writeheader()

    genotype_counter = itertools.count(0)
    genotypes = create_population(population_size=args.population_size,
                                  counter=genotype_counter)

    for gen in srange(args.generations):
        genotypes.vs.select(lambda v: v['first_seen'] is None)['first_seen'] = gen
        genotypes.vs.select(lambda v: v['abundance'] > 0)['last_seen'] = gen
        genotypes.vs.select(lambda v: v['total_abundance'] > 0)['lineage_last_seen'] = gen

        if not args.quiet:
            sys.stdout.write("\r")
            pct = float(gen) / args.generations
            sys.stdout.write("[%-40s] %d%%" % ('='* int(40 * pct), (pct * 100)))
            sys.stdout.flush()

        for g in genotypes.vs.select(lambda v: v['total_abundance'] > 0):
            outfile.writerow({'Generation': gen,
                              'Genotype': g['name'],
                              'Parent': g['parent'],
                              'FirstSeen': g['first_seen'],
                              'LastSeen': g['last_seen'],
                              'LineageLastSeen': g['lineage_last_seen'],
                              'Depth': g['depth'],
                              'Fitness': g['fitness'],
                              'FitnessEffects': sum(g['fitness_effects']),
                              'FitnessDiff': g['fitness_diff'],
                              'Abundance': g['abundance'],
                              'TotalAbundance': g['total_abundance'],
                              'Frequency': g['frequency'],
                              'MaxFrequency': g['max_frequency'],
                              'TotalFrequency': g['total_frequency'],
                              'FixationTime': g['fixation_time']})

        reproduce(population=genotypes, population_size=args.population_size)
        #mutate_multiples(population=genotypes,
        #                 mutation_rate=args.mutation_rate,
        #                 genotype_counter=genotype_counter)
        mutate_bdc(population=genotypes, mutation_rate=args.mutation_rate,
                   genotype_counter=genotype_counter)

        prune_frequency(population=genotypes, min_frequency=args.prune_freq)
        get_total_abundances(genotype=genotypes.vs[0])

        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None):
            v['abundances'].append(int(v['abundance']))

        genotypes['generations'] += 1
        genotypes['population_size'] = int(sum(genotypes.vs['abundance']))
        genotypes.vs['frequency'] = genotypes.vs['abundance'] / sum(genotypes.vs['abundance'])
        genotypes.vs['total_frequency'] = genotypes.vs['total_abundance'] / sum(genotypes.vs['abundance'])


        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']

        genotypes.vs.select(lambda v: v['fixation_time'] == -1 and (float(v['total_abundance']) / genotypes['population_size']) >= args.fixation_freq)['fixation_time'] = gen

        # For writing the tree at every cycle
        #genotypes.write_gml("TREES/genotypes-{0:06d}.gml".format(gen))

    graph_write_csv(genotypes, os.path.join(args.data_dir, "tree-end.csv"))
    graph_write_json(genotypes, os.path.join(args.data_dir, "tree-end.json"),
                     sort_keys=True)
    genotypes.write_gml(os.path.join(args.data_dir, "tree-end.gml"))

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    run_simulation()

