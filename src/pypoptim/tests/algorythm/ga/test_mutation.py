import pytest

from ....algorythm.ga.mutation import cauchy_inverse_cdf, cauchy_mutation, cauchy_mutation_population, reflection
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


def test_cauchy_mutation_population(population):

    bounds = [[-5, 5], [3, 13]]
    n_organisms = 0
    p = population(n_organisms, bounds=bounds)
    new_p = cauchy_mutation_population(population=p, bounds=bounds, gamma=1, mutation_rate=1)
    assert new_p==[]

    n_organisms = 2
    p = population(n_organisms, bounds=bounds)
    with pytest.raises(ValueError):
        cauchy_mutation_population(population=p, bounds=bounds * 2, gamma=1, mutation_rate=0)

    for g, mut in zip([0, 1], [1, 0]):
        new_p = cauchy_mutation_population(population=p, bounds=bounds, gamma=g, mutation_rate=mut)
        for organism, new_organism in zip(p, new_p):
            new_organism.update()
            assert organism.x.all() == new_organism.x.all()

    seed = 888
    new_p = [cauchy_mutation_population(population=p,
                                        bounds=bounds,
                                        gamma=1,
                                        mutation_rate=1,
                                        rng= np.random.default_rng(seed)) for _ in range(2)]
    for org_1, org_2 in zip(new_p[0], new_p[1]):
        org_1.update()
        org_2.update()
        assert org_1.x.all() == org_2.x.all()

    assert new_p[0] == new_p[1]

def test_reflection():
    ub = np.array([-1, -1, -1, -2])
    lb = np.array([1, 1, 1, 1])
    genes = np.array([0.5, 0.5, 0.5, 0])
    shifts = np.array([0.5, 0.8, 1.3, -5.5])
    assert reflection(ub=ub, lb=lb, genes=genes, shifts=shifts).all() == np.array([1, 0.3, 0.2, 0.5]).all()



