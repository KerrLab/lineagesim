import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import exponential as nexp
from scipy.special import gammainc as gamma
from six.moves import range as srange


def mutate_multiples(population, mutation_rate, genotype_counter):
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
        num_k_mutants = np.subtract(r, new_r) #number of k-mutants
        for parent_id in np.nonzero(num_k_mutants)[0]:
            for _mutant in srange(num_k_mutants[parent_id]):
                mu_effect = nexp(scale=0.01, size=k)
                population.add_vertex(name=next(genotype_counter),
                                      parent=int(population.vs[parent_id]['name']),
                                      abundance=1,
                                      abundances=[1],
                                      total_abundance=1,
                                      depth=population.vs[parent_id]['depth'] + 1,
                                      fitness=population.vs[parent_id]['fitness'] + np.prod(1 + mu_effect),
                                      fitness_effects=mu_effect.tolist(),
                                      fitness_diff=np.prod(1 + mu_effect),
                                      frequency=1 / population['population_size'],
                                      max_frequency=0,
                                      fixation_time=-1)
                population.add_edge(source=parent_id,
                                    target=population.vcount() - 1,
                                    fitness_effect=np.prod(1 + mu_effect))

        r = new_r
        k += 1
    return population

