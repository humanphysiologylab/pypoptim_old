# import numpy as np

from ..algorythm.ga.selection import tournament_selection


def test_tournament_selection():
    p = [1, 0]
    assert tournament_selection(p) is p[1]
    p.append(-1)
    assert tournament_selection(p, k=3) is p[-1]

#
# def test_do_step():
#
#     def create_population(n_organisms=3):
#
#         population = []
#         bounds = []
#         gammas = []
#
#         for i in range(n_organisms):
#             organism = dict(genes=np.full(n_organisms, i),
#                             state=str(i),
#                             fitness=i)
#             population.append(organism)
#
#             bounds.append([0 - 1, n_organisms])
#             gammas.append(1)
#
#         bounds = np.array(bounds)
#         gammas = np.array(gammas)
#
#         return population, bounds, gammas
#
#     kw = dict(crossover_rate=1.,
#               mutation_rate=1.,
#               gamma=1.)
#
#     n_organisms = 3
#
#     population, bounds, gammas = create_population(n_organisms)
#     new_population = ga.do_step(population=population,
#                                 new_size=n_organisms,
#                                 elite_size=n_organisms,
#                                 bounds=bounds,
#                                 **kw)
#     assert np.all(organism in new_population for organism in population)
#
#     population, bounds, gammas = create_population(n_organisms)
#     new_population = ga.do_step(population=population,
#                                 new_size=0,
#                                 elite_size=n_organisms,
#                                 bounds=bounds,
#                                 **kw)
#     assert len(new_population) == 0
#
#     population, bounds, gammas = create_population(n_organisms)
#     new_population = ga.do_step(population=population,
#                                 new_size=101,
#                                 elite_size=0,
#                                 bounds=bounds,
#                                 **kw)
#     assert all(organism['state'] != '0' for organism in new_population)
