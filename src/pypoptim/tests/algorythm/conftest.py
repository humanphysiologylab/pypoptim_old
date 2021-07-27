import pytest
import numpy as np
from ...algorythm.solution import Solution
from ...algorythm.ga import GA


class SquareSolution(Solution):
    def update(self):
        self._y = np.sum(self.x ** 2)

    def is_valid(self):
        return self.is_updated()


class MaxAbsSolution(Solution):
    def update(self):
        self._y = np.max(np.abs(self.x))

    def is_valid(self):
        return self.is_updated()


@pytest.fixture()
def square_solution():
    return SquareSolution


@pytest.fixture()
def maxabs_solution():
    return MaxAbsSolution


@pytest.fixture()
def population():
    def _population(n, bounds, **kw):
        ga_optim = GA(SolutionSubclass=SquareSolution, bounds=bounds, **kw)
        return ga_optim.generate_population(n)
    return _population
