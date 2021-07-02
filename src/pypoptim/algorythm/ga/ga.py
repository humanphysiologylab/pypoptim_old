import copy
import numpy as np

from pypoptim.helpers import random_value_from_bounds
import random
from .crossover import sbx_crossover
from .mutation import cauchy_mutation_population
from ..solution import Solution


class GA:

    def __init__(self, SolutionSubclass, bounds, gammas=None,
                 mutation_rate=1., crossover_rate=1., selection_force=2,
                 keys_data_transmit=None):

        if not issubclass(SolutionSubclass, Solution):
            raise TypeError
        self._SolutionSubclass = SolutionSubclass

        if not isinstance(bounds, np.ndarray):
            raise TypeError
        if bounds.ndim != 2 or bounds.shape[0] == 0 or bounds.shape[1] != 2:
            raise ValueError
        if np.any(bounds[:, 0] >= bounds[:, 1]):
            raise ValueError
        self._bounds = bounds
        self._n_genes = len(bounds)

        self.__gamma_default = 1
        if gammas is None:
            self._gammas = np.full(len(bounds), self.__gamma_default)
        else:
            if not isinstance(gammas, np.ndarray):
                raise TypeError
            if len(gammas) != self._n_genes:
                raise ValueError
            if np.any(gammas <= 0):
                raise ValueError
            self._gammas = gammas

        if not (0. <= mutation_rate <= 1):
            raise ValueError
        self._mutation_rate = mutation_rate

        if not (0. <= crossover_rate <= 1):
            raise ValueError
        self._crossover_rate = crossover_rate

        if not isinstance(selection_force, int):
            raise TypeError
        if selection_force < 2:
            raise ValueError
        self._selection_force = selection_force

        if keys_data_transmit is not None:
            if not isinstance(keys_data_transmit, (list, tuple)):
                raise ValueError
            self._keys_data_transmit = keys_data_transmit
        else:
            self._keys_data_transmit = []

        self._population = []

    def __repr__(self):
        s =  f'{self.__class__.__name__}:\n'
        s += f'    bounds = {self._bounds}\n'
        s += f'    gammas = {self._gammas}\n'
        s += f'    mutation_rate   = {self._mutation_rate}\n'
        s += f'    crossover_rate  = {self._crossover_rate}\n'
        s += f'    selection_force = {self._selection_force}\n'
        s += f'    keys_data_transmit = {self._keys_data_transmit}\n'
        return s

    def __str__(self):
        return self.__repr__()

    def generate_solution(self) -> Solution:
        genes = [random_value_from_bounds(self._bounds[i]) for i in range(self._n_genes)]
        sol = self._SolutionSubclass()
        sol.x = np.array(genes)
        return sol

    def generate_population(self, n_solutions: int) -> None:
        self._population = [self.generate_solution() for _ in range(n_solutions)]

    def _transmit_solution_data(self, sol_parent: Solution, sol_child: Solution):
        for key in self._keys_data_transmit:
            if key not in sol_parent:
                raise KeyError
            sol_child[key] = sol_parent[key]

    def _selection(self) -> Solution:  # tournament selection
        return min(random.sample(self._population, self._selection_force))

    def _crossover(self, genes1, genes2) -> tuple:
        return sbx_crossover(genes1, genes2, bounds=self._bounds)

    def get_mutants(self, size=1):

        if not isinstance(size, int):
            raise TypeError
        if size < 0:
            raise ValueError

        new_population = []

        while len(new_population) < size:
            if not len(self._population):
                raise ValueError
            parent1, parent2 = self._population[0], self._population[0]
            while parent1 is parent2:
                parent1 = self._selection()
                parent2 = self._selection()

            if np.random.random() <= self._crossover_rate:
                offspring_genes = self._crossover(parent1.x, parent2.x)
                parent_data_transmitter = min(parent1, parent2)
                for genes in offspring_genes:
                    child = self._SolutionSubclass()
                    child.x = genes
                    self._transmit_solution_data(parent_data_transmitter, child)  # child['state'] = parent1['state']
                    new_population.append(child)
            else:  # no crossover
                child1 = copy.deepcopy(parent1)
                child2 = copy.deepcopy(parent2)
                new_population += [child1, child2]

        new_population = new_population[:size]  # TODO: sbx_crossover creates pairs so this is for odd size of the population

        cauchy_mutation_population(new_population,
                                   bounds=self._bounds,
                                   gamma=self.__gamma_default,
                                   mutation_rate=self._mutation_rate,
                                   inplace=True)
        return new_population

    def get_elites(self, size=1):
        if not isinstance(size, int):
            raise TypeError
        if not (0 <= size < len(self._population)):
            raise ValueError
        return sorted(self._population)[:size]

    @property
    def population(self):
        return copy.deepcopy(self._population)

    @population.setter
    def population(self, population):
        self._population = copy.deepcopy(population)

    def update_population(self) -> None:
        for sol in self._population:
            sol.update()


    def run(self, n_solutions, n_epochs=1, n_elites=0):

        if not isinstance(n_solutions, int):
            raise TypeError
        if n_solutions < 0:
            raise ValueError

        if not isinstance(n_epochs, int):
            raise TypeError
        if n_epochs < 0:
            raise ValueError

        if not isinstance(n_elites, int):
            raise TypeError
        if not (0 <= n_elites <= n_solutions):
            raise ValueError

        if not self._population:
            self.generate_population(n_solutions)
        for i in range(n_epochs):
            self.update_population()
            elites = self.get_elites(n_elites)
            mutants = self.get_mutants(n_solutions - n_elites)
            self._population = elites + mutants
