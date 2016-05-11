import igraph
import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import multinomial as nmultinom
from numpy.random import normal as nnormal
from six.moves import range as srange

import csv
import itertools

POPSIZE = int(1e6)
NUM_CYCLES = 1000
MUTATION_RATE = 1e-6
OUTFILENAME = "results.csv"
THRESH_FREQ = 0.01


class Population(object):
    """A population object stores the different genotypes present in the population and their abundances

    Each genotype is represented as a node on a tree, which contains additional
    information about that genotype. When a mutation occurs, a daughter node
    is created.
    """
    def __init__(self, population_size, base_fitness = 1.0):
        """Create a population"""
        self.genotype_counter = itertools.count(0)
        self.conf_popsize = population_size
        self.generations = 0
        self.genotypes = igraph.Graph(directed = True,
                                      graph_attrs = {'population_size': population_size,
                                                     'generations': 0},
                                      vertex_attrs = {'first_seen': None,
                                                      'last_seen': None,
                                                      'depth': None,
                                                      'abundance': 0,
                                                      'abundances': [0],
                                                      'frequency': 0,
                                                      'max_frequency': 0,
                                                      'fitness': 0,
                                                      'fitness_diff': [0]})

        self.genotypes.add_vertex(name = next(self.genotype_counter), depth = 0,
                                  abundance = population_size,
                                  abundances = [population_size],
                                  frequency = 1.0, max_frequency = 1.0,
                                  fitness = base_fitness,
                                  fitness_diff = [0])

    def set_seen_times(self, time):
        """Update each genotype's first- and last seen times"""
        self.genotypes.vs.select(lambda v: v['first_seen'] is None)['first_seen'] = time
        self.genotypes.vs.select(lambda v: v['abundance'] > 0)['last_seen'] = time

    def set_frequencies(self):
        """Update each genotype's frequency and maximum frequency"""
        self.genotypes.vs['frequency'] = self.genotypes.vs['abundance'] / sum(self.genotypes.vs['abundance'])

        for v in self.genotypes.vs.select(lambda v: v['abundance'] > 0 and v['first_seen'] is not None and v['frequency'] > v['max_frequency']):
            v['max_frequency'] = v['frequency']

    def reproduce(self):
        """Fill the population with genotypes proportional to their fitness"""
        fitnesses = np.array(self.genotypes.vs['fitness'])
        abundances = np.array(self.genotypes.vs['abundance'])
        ab_fit = fitnesses * abundances

        self.genotypes.vs['abundance'] = nmultinom(n = self.conf_popsize,
                                                   pvals = ab_fit / ab_fit.sum(),
                                                   size = 1)[0]

    def mutate(self, mutation_rate):
        """Mutate individuals in the population, creating new genotypes"""
        assert(mutation_rate >= 0 and mutation_rate <= 1)

        num_mutants = nbinom(n = self.genotypes.vs['abundance'],
                             p = mutation_rate)
        self.genotypes.vs['abundance'] = self.genotypes.vs['abundance'] - num_mutants

        for parent_id in np.nonzero(num_mutants)[0]:
            for _mutant in srange(num_mutants[parent_id]):

                mu_effect = nnormal(loc = 0.0, scale = 0.1)

                self.genotypes.add_vertex(name = next(self.genotype_counter),
                                          abundance = 1, abundances = [1],
                                          depth = self.genotypes.vs[parent_id]['depth'] + 1,
                                          fitness = self.genotypes.vs[parent_id]['fitness'] + mu_effect,
                                          fitness_diff = [mu_effect],
                                          frequency = 0, max_frequency = 0)
                self.genotypes.add_edge(source = parent_id,
                                        target = self.genotypes.vcount() - 1,
                                        color = int(mu_effect >= 0))


    def dilute(self, dilution_prob):
        """Probabilistically thin the population"""
        assert(dilution_prob >= 0 and dilution_prob <= 1)
        self.genotypes.vs['abundance'] = nbinom(n = self.genotypes.vs['abundance'],
                                                p = dilution_prob)

    def prune(self, min_frequency):
        """Prune extinct leaf nodes that did not reach a given frequency"""
        self.genotypes.vs.select(lambda v: v.outdegree() == 0 and v['abundance'] == 0 and v['first_seen'] is not None and v['max_frequency'] < min_frequency).delete()

    def evolve(self, time, mutation_rate, prune_min_freq):
        self.set_seen_times(time = time)
        self.reproduce()
        self.mutate(mutation_rate = mutation_rate)
        self.set_frequencies()
        self.prune(min_frequency = prune_min_freq)
        self.generations += 1

    def write_gml(self, filename, **kwargs):
        """Write a GML file containing the current state of the population"""
        self.genotypes.write_gml(f = filename, **kwargs)

    def write_json(self, filename):
        """Write a JSON file containing the current state of the population"""
        pass

    def write_csvdata(self, outfile, time):
        """Write information about each genotype to a CSV file"""
        for extant in self.genotypes.vs.select(lambda v: v['abundance'] > 0):
            outfile.writerow({'Generation': time,
                              'Genotype': extant.index,
                              'Depth': extant['depth'],
                              'Fitness': extant['fitness'],
                              'Abundance': extant['abundance'],
                              'Frequency': extant['frequency']})

# ---------------

def run_simulation(num_generations, outfile):
    outfile = csv.DictWriter(open(outfile, 'w'),
                             fieldnames = ['Generation', 'Genotype', 'Depth', 'Fitness', 'Abundance', 'Frequency'])
    outfile.writeheader()

    p = Population(population_size = POPSIZE)

    for gen in srange(num_generations):
        print("Gen {g}".format(g = gen))
        p.evolve(time = gen, mutation_rate = MUTATION_RATE,
                 prune_min_freq = THRESH_FREQ)
        p.write_csvdata(outfile = outfile, time = gen)

    p.write_gml("tree-end.gml")


if __name__ == "__main__":
    run_simulation(num_generations = NUM_CYCLES, outfile = OUTFILENAME)
