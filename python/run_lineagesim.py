#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import datetime
import csv
import itertools
import os
import sys
from warnings import warn

import numpy as np

from six.moves import range as srange

import lineagesim as ls


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
                        default=1.0, type=check_positive_01,
                        help='Threshold frequency for classification as fixed (default: 1)')
    parser.add_argument('--generations', '-G', metavar='G', default=1000,
                        type=check_positive_int,
                        help='Number of generations to simulate (default: 1000)')
    parser.add_argument('--mutation_rate', '-m', metavar='μ', default=1e-6,
                        type=check_positive_01,
                        help='Mutation rate (default: 1e-6)')
    parser.add_argument('--mutation_scale', '-s', metavar='s', default=0.01,
                        type=check_positive_01,
                        help='Mutation scale (default: 0.01)')
    parser.add_argument('--population_size', '-N', metavar='N',
                        default=int(1e6), type=float,
                        help='Size of the population (default: 1,000,000)')
    parser.add_argument('--prune_freq', '-p', metavar='F',
                        type=check_positive_01,
                        help='Threshold frequency for pruning (default: (N * μ) / N)')
    parser.add_argument('--seed', '-S', metavar='S', help='Set the '\
                        'pseudorandom number generator seed',
                        type=check_nonnegative_int)
    parser.add_argument('--quiet', '-q', action='store_true', default=False,
                        help='Suppress output messages')
    parser.add_argument('--version', action='version', version=ls.__version__)

    return parser.parse_args()


def run_simulation(args=parse_arguments()):
    """Run the simulation"""

    # Make sure population size is an integer. Doing the conversion here so
    # that large numbers like 1e7, which are floats, can be given as args.
    if int(args.population_size) != args.population_size:
        msg = "Population size must be integer. Using {N}.".format(N=int(args.population_size))
        warn(msg)
    args.population_size = int(args.population_size)

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

    outfile = csv.DictWriter(open(os.path.join(args.data_dir, "lineages.csv"),
                                  'w'),
                             fieldnames=['Generation',
                                         'Genotype',
                                         'Parent',
                                         'FirstSeen',
                                         'LastSeen',
                                         'LineageLastSeen',
                                         'Depth',
                                         'Fitness',
                                         'FitnessDiff',
                                         'NumMutations',
                                         'Abundance',
                                         'LineageAbundance',
                                         'Frequency',
                                         'MaxFrequency',
                                         'LineageFrequency',
                                         'FixationTime'])
    outfile.writeheader()

    genotype_counter = itertools.count(0)

    genotypes = ls.create_population(population_size=int(args.population_size),
                                     counter=genotype_counter)

    for gen in srange(args.generations):
        genotypes.vs.select(lambda v: v['first_seen'] is None)['first_seen'] = gen
        genotypes.vs.select(lambda v: v['abundance'] > 0)['last_seen'] = gen
        genotypes.vs.select(lambda v: v['lineage_abundance'] > 0)['lineage_last_seen'] = gen

        if not args.quiet:
            sys.stdout.write("\r")
            pct = float(gen) / args.generations
            sys.stdout.write("[%-40s] %d%%" % ('='* int(40 * pct), (pct * 100)))
            sys.stdout.flush()

        for g in genotypes.vs.select(lambda v: v['lineage_abundance'] > 0):
            outfile.writerow({'Generation': gen,
                              'Genotype': g['name'],
                              'Parent': g['parent'],
                              'FirstSeen': g['first_seen'],
                              'LastSeen': g['last_seen'],
                              'LineageLastSeen': g['lineage_last_seen'],
                              'Depth': g['depth'],
                              'Fitness': g['fitness'],
                              'FitnessDiff': g['fitness_diff'],
                              'NumMutations': len(g['fitness_effects']),
                              'Abundance': g['abundance'],
                              'LineageAbundance': g['lineage_abundance'],
                              'Frequency': g['frequency'],
                              'MaxFrequency': g['max_frequency'],
                              'LineageFrequency': g['lineage_frequency'],
                              'FixationTime': g['fixation_time']})

        ls.reproduce(population=genotypes,
                     population_size=args.population_size)
        ls.mutate(population=genotypes,
                  mutation_rate=args.mutation_rate,
                  mutation_scale=args.mutation_scale,
                  genotype_counter=genotype_counter)

        ls.prune_frequency(population=genotypes, min_frequency=args.prune_freq)
        ls.get_lineage_abundances(genotype=genotypes.vs[0])

        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None):
            v['abundances'].append(int(v['abundance']))

        genotypes['generations'] += 1
        genotypes['population_size'] = int(sum(genotypes.vs['abundance']))
        genotypes.vs['frequency'] = genotypes.vs['abundance'] / sum(genotypes.vs['abundance'])
        genotypes.vs['lineage_frequency'] = genotypes.vs['lineage_abundance'] / sum(genotypes.vs['abundance'])

        for v in genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']

        genotypes.vs.select(lambda v: v['fixation_time'] == -1 and (float(v['lineage_abundance']) / genotypes['population_size']) >= args.fixation_freq)['fixation_time'] = gen

        # For writing the tree at every cycle
        #genotypes.write_gml("TREES/genotypes-{0:06d}.gml".format(gen))

    ls.output.graph_write_csv(genotypes,
                              os.path.join(args.data_dir, "lineages-end.csv"))
    ls.output.graph_write_json(genotypes,
                               os.path.join(args.data_dir, "tree-end.json"),
                               sort_keys=True)
    genotypes.write_gml(os.path.join(args.data_dir, "tree-end.gml"))



if __name__ == "__main__":
    run_simulation()

