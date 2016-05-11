import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import exponential as nexp
from scipy.special import gammainc as gamma
from six.moves import range as srange

def mutate_multiples(p, mutation_rate, counter):
    """Mutate individuals in the population"""

    assert(mutation_rate >= 0)

    num_mutants = nbinom(n = p.vs['abundance'], p = 1 - np.exp(-mutation_rate))
    p.vs['abundance'] = p.vs['abundance'] - num_mutants

    r = num_mutants
    k = 1
    while any(r > 0):
        #r is the number of individuals with at least k+1 mutations
        #rate inside binomial is the prob of an individual having at least k+1 mutations given that it has at least k mutations
        new_r = nbinom(n = r, p = gamma(k+1, mutation_rate)/gamma(k, mutation_rate), size = len(r))
        num_k_mutants = np.subtract(r,new_r) #number of k-mutants
        for parent_id in np.nonzero(num_k_mutants)[0]:
            for _mutant in srange(num_k_mutants[parent_id]):
                mu_effect = nexp(scale = 0.01, size = k)
                p.add_vertex(name = next(counter), abundance = 1,
                             abundances = [1],
                             depth = p.vs[parent_id]['depth'] + 1,
                             fitness = p.vs[parent_id]['fitness'] + np.prod(1 + mu_effect),
                             fitness_diff = [mu_effect], frequency = 0,
                             max_frequency = 0)
                p.add_edge(source = parent_id, target = p.vcount() - 1,
                           color = int(np.prod(1 + mu_effect) >= 0))

        r = new_r
        k += 1
    return p
