from ....algorythm.ga.mutation import cauchy_inverse_cdf, cauchy_mutation, cauchy_mutation_population
import numpy as np


def test_cauchy_inverse_cdf():

    gamma = 0
    assert cauchy_inverse_cdf(gamma=gamma, rng=None) == 0

    gamma = 1
    seed = 888
    x = [cauchy_inverse_cdf(gamma=gamma, rng=np.random.default_rng(seed)) for _ in range(2)]
    assert x[0] == x[1]


def test_cauchy_mutation():
    genes = [-1, 1]
    gamma = 0
    assert cauchy_mutation(genes=genes, gamma=gamma) == genes

    seed = 888
    x = [cauchy_mutation(genes=genes, rng=np.random.default_rng(seed)) for _ in range(2)]
    assert x[0] == x[1]




def test_cauchy_mutation_population():
    assert 0
