import pytest
import numpy as np
from ....algorythm.ga.ga import GA


@pytest.fixture()
def ga_optim_fabric(square_solution):
    def _ga_optim_fabric(**kw):
        bounds = np.asfarray([[-1, 1],
                              [2, 4]])
        rng = np.random.default_rng(42)
        return GA(SolutionSubclass=square_solution,
                  bounds=bounds,
                  rng=rng,
                  **kw)
    return _ga_optim_fabric


@pytest.fixture()
def ga_optim_default(ga_optim_fabric):
    return ga_optim_fabric()


@pytest.fixture()
def ga_optim_with_data(ga_optim_fabric):
    keys_data_transmit = ['state']
    ga_optim_fabric(keys_data_transmit=keys_data_transmit)

@pytest.fixture()
def ga_optim_for_is_valid(ga_optim_fabric):
    ga_optim = ga_optim_fabric()
    class SS(ga_optim._SolutionSubclass):
        def is_valid(self):
            return np.all(self.x > 0)
    ga_optim._SolutionSubclass = SS
    return ga_optim


class TestGA:

    def test_init(self, square_solution):
        assert 0

    def test_is_inside_bounds(self, ga_optim_default):
        ga_optim = ga_optim_default
        bounds = ga_optim.bounds

        sol = ga_optim_default.generate_solution()
        sol.x = [0, 3]
        assert ga_optim._is_solution_inside_bounds(sol, bounds)
        assert ga_optim._is_solution_inside_bounds(sol)

        sol.x = [0, 0]
        assert not ga_optim._is_solution_inside_bounds(sol, bounds)
        assert not ga_optim._is_solution_inside_bounds(sol)

    def test_generate_solution(self, ga_optim_fabric):

        ga_optim_1 = ga_optim_fabric()
        sol_1 = ga_optim_1.generate_solution()

        ga_optim_2 = ga_optim_fabric()
        sol_2 = ga_optim_2.generate_solution()
        sol_2_another = ga_optim_2.generate_solution()

        assert np.all(sol_1.x == sol_2.x)
        assert np.all(sol_2.x != sol_2_another.x)  # super likely
        assert not sol_1.is_updated()
        assert not sol_1.is_valid()

        for sol in sol_1, sol_2, sol_2_another:
            assert ga_optim_1._is_solution_inside_bounds(sol)

    def test_generate_population(self, ga_optim_default):
        n = 42
        population = ga_optim_default.generate_population(n)
        assert len(population) == n
        assert all(ga_optim_default._is_solution_inside_bounds(sol) for sol in population)

    def test_update_population(self, ga_optim_default):
        n = 42
        population = ga_optim_default.generate_population(n)
        assert all(not sol.is_updated() for sol in population)

        ga_optim_default.update_population(population)
        assert all(sol.is_updated() for sol in population)

    def test_filter_population(self, ga_optim_for_is_valid):
        n = 42
        population = ga_optim_for_is_valid.generate_population(n)

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population_filtered) == 0  # no one is updated

        bounds = ga_optim_for_is_valid.bounds
        margin = 0.1
        assert np.all(bounds[:, 1] > margin)
        for sol in population:
            sol.x = bounds[:, 1] - margin
            sol.update()
            assert ga_optim_for_is_valid._is_solution_inside_bounds(sol)

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population) == len(population_filtered)
        assert population == population_filtered

        population[0].x = population[0].x  # makes this solution not updated
        assert not population[0].is_updated()
        population[1]._x = np.full_like(len(bounds), -0.1)  # makes this solution invalid, yet updated
        assert not population[1].is_valid() and population[1].is_updated()
        population[2]._x = bounds[:, 1] + 1  # makes this solution outside bounds, yet updated and valid
        assert not ga_optim_for_is_valid._is_solution_inside_bounds(population[2])\
               and population[2].is_updated() and population[2].is_valid()

        population_filtered = ga_optim_for_is_valid.filter_population(population)
        assert len(population) == len(population_filtered) + 3
        assert population[3:] == population_filtered

    def test_transmit_solution_data(self):
        assert 0

    def test_get_mutants(self):
        assert 0

    def test_get_elites(self, ):
        assert 0

    def test_run(self):
        assert 0
