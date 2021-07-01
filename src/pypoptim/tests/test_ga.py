import numpy as np

from ..algorythm import ga


def test_transform_genes():

    genes  = np.array( [1.,      1.0,    -1.,   100.0])

    bounds = np.array([[0.,      0.1,    -3.,    10.0],
                       [2.,     10.0,     0.,  1000.0]]).T

    gammas = np.array( [1.,      1.0,     2.,     0.5])

    mask_multipliers = np.array([False, True, False, True])

    genes_transformed, bounds_transformed = ga.transform_genes_bounds(genes, bounds, gammas, mask_multipliers)

    assert np.allclose(genes_transformed, np.array([0.5, 0.5, 0.5, 1]))

    assert np.allclose(bounds_transformed, np.array([[0.  , 1.  ],
                                                     [0.  , 1.  ],
                                                     [0.  , 0.75],
                                                     [0.  , 2.  ]]))

    genes_back = ga.transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds, mask_multipliers)

    assert np.allclose(genes, genes_back)


def test_tournament_selection():

    p = [dict(dummy=0, fitness=1),
         dict(dummy=1, fitness=0)]

    assert ga.tournament_selection(p) is p[0]

    p.append(dict(foo=42, bar='bar', fitness=100))

    assert ga.tournament_selection(p, k=3) is p[-1]


def test_do_step():

    def create_population(n_organisms=3):

        population = []
        bounds = []
        gammas = []

        for i in range(n_organisms):
            organism = dict(genes=np.full(n_organisms, i),
                            state=str(i),
                            fitness=i)
            population.append(organism)

            bounds.append([0 - 1, n_organisms])
            gammas.append(1)

        bounds = np.array(bounds)
        gammas = np.array(gammas)

        return population, bounds, gammas

    kw = dict(crossover_rate=1.,
              mutation_rate=1.,
              gamma=1.)

    n_organisms = 3

    population, bounds, gammas = create_population(n_organisms)
    new_population = ga.do_step(population=population,
                                new_size=n_organisms,
                                elite_size=n_organisms,
                                bounds=bounds,
                                **kw)
    assert np.all(organism in new_population for organism in population)

    population, bounds, gammas = create_population(n_organisms)
    new_population = ga.do_step(population=population,
                                new_size=0,
                                elite_size=n_organisms,
                                bounds=bounds,
                                **kw)
    assert len(new_population) == 0

    population, bounds, gammas = create_population(n_organisms)
    new_population = ga.do_step(population=population,
                                new_size=101,
                                elite_size=0,
                                bounds=bounds,
                                **kw)
    assert all(organism['state'] != '0' for organism in new_population)
