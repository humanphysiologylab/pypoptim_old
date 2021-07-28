import copy

import numpy as np

from pypoptim.helpers import random_value_from_bounds,\
                             transform_genes_bounds,\
                             transform_genes_bounds_back, \
                             is_values_inside_bounds
from .selection import tournament_selection
from .crossover import sbx_crossover
from .mutation import cauchy_mutation_population
from ..solution import Solution

import logging
logger = logging.getLogger(__name__)

class GA:

    def __init__(self, SolutionSubclass, bounds, gammas=None, mask_log10_scale=None,
                 mutation_rate=1., crossover_rate=1., selection_force=2,
                 keys_data_transmit=None, rng=None):

        if issubclass(SolutionSubclass, Solution) and not isinstance(SolutionSubclass, Solution):
            self._SolutionSubclass = SolutionSubclass
        else:
            raise TypeError

        bounds = np.asfarray(bounds)
        if bounds.ndim != 2 or bounds.shape[0] == 0 or bounds.shape[1] != 2:
            raise ValueError
        if np.any(bounds[:, 0] >= bounds[:, 1]):
            raise ValueError
        self._bounds = bounds
        self._n_genes = len(bounds)

        self.__gamma_default = 1
        if gammas is None:
            self._gammas = np.full(self._n_genes, self.__gamma_default)
        else:
            gammas = np.asfarray(gammas)
            if len(gammas) != self._n_genes:
                raise ValueError
            if np.any(gammas <= 0):
                raise ValueError
            self._gammas = gammas

        if mask_log10_scale is None:
            self._mask_log10_scale = np.full(self._n_genes, False)
        else:
            mask_log10_scale = np.asarray(mask_log10_scale)
            if len(mask_log10_scale) != self._n_genes:
                raise ValueError
            if np.any(self._bounds[mask_log10_scale.astype(bool), 0] <= 0):
                raise ValueError
            self._mask_log10_scale = mask_log10_scale

        _, self._bounds_transformed = transform_genes_bounds(self._bounds[:, 0],
                                                             self._bounds,
                                                             self._gammas,
                                                             self._mask_log10_scale)

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
                raise TypeError
            self._keys_data_transmit = keys_data_transmit
        else:
            self._keys_data_transmit = []

        if rng is not None:
            if not isinstance(rng, np.random._generator.Generator):
                raise TypeError
            self._rng = rng
        else:
            self._rng = np.random.default_rng()

    @property
    def bounds(self):
        return self._bounds.copy()

    def __repr__(self):

        hstack = np.hstack([self._bounds,
                            self._gammas.reshape(-1, 1),
                            self._mask_log10_scale.reshape(-1, 1)])
        s =  f'{self.__class__.__name__}:\n'
        s += f'[bounds_lower bounds_upper gammas mask_log_10_scale]\n'
        s += str(hstack) + "\n"
        s += f'mutation_rate: {self._mutation_rate}\n'
        s += f'crossover_rate: {self._crossover_rate}\n'
        s += f'selection_force: {self._selection_force}\n'
        s += f'keys_data_transmit: {self._keys_data_transmit}\n'
        return s

    def __str__(self):
        return self.__repr__()

    def _transform_genes(self, genes):
        genes_transformed, bounds_transformed = transform_genes_bounds(genes,
                                                                       self._bounds,
                                                                       self._gammas,
                                                                       self._mask_log10_scale)
        assert(np.allclose(bounds_transformed, self._bounds_transformed))
        return genes_transformed

    def _transform_genes_back(self, genes_transformed):
        genes = transform_genes_bounds_back(genes_transformed,
                                            self._bounds_transformed,
                                            self._bounds,
                                            self._mask_log10_scale)
        return genes

    def _transform_solution(self, sol):
        sol_transformed = self._SolutionSubclass(self._transform_genes(sol.x), **sol.data)
        sol_transformed._y = sol.y
        return sol_transformed

    def _transform_solution_back(self, sol_transformed):
        genes = self._transform_genes_back(sol_transformed.x)
        sol = self._SolutionSubclass(genes, **sol_transformed.data)
        sol._y = sol_transformed.y
        return sol

    def generate_solution(self) -> Solution:
        genes_transformed = [random_value_from_bounds(self._bounds_transformed[i], rng=self._rng)
                             for i in range(self._n_genes)]
        genes_transformed = np.array(genes_transformed)
        genes = self._transform_genes_back(genes_transformed)
        sol = self._SolutionSubclass(genes)
        return sol

    def generate_population(self, n_solutions: int) -> list:
        return [self.generate_solution() for _ in range(n_solutions)]

    def _transmit_solution_data(self, sol_parent: Solution, sol_child: Solution):
        for key in self._keys_data_transmit:
            sol_child[key] = copy.deepcopy(sol_parent[key])

    def _selection(self, population) -> Solution:  # tournament selection
        return tournament_selection(np.array(population, dtype=object),
                                    self._selection_force, rng=self._rng)

    def _crossover(self, genes1, genes2) -> tuple:
        return sbx_crossover(genes1, genes2,
                             bounds=self._bounds_transformed,
                             cross_rate=self._crossover_rate,
                             rng=self._rng)

    def get_mutants(self, population, size=1):

        if len(population) < 3:
            raise ValueError("Insufficient number of solution for the crossover")
        if not isinstance(size, int):
            raise TypeError
        if size < 0:
            raise ValueError

        new_population = []

        while len(new_population) < size:

            parent1, parent2 = population[0], population[0]
            while parent1 is parent2:
                parent1 = self._selection(population)
                parent2 = self._selection(population)

            parent1_transformed, parent2_transformed = [self._transform_solution(p) for p in (parent1, parent2)]

            if self._rng.random() <= self._crossover_rate:
                offspring_genes = self._crossover(parent1_transformed.x,
                                                  parent2_transformed.x)
                parent_data_transmitter = min(parent1_transformed, parent2_transformed)
                for genes in offspring_genes:
                    child = self._SolutionSubclass(genes)
                    self._transmit_solution_data(parent_data_transmitter, child)
                    new_population.append(child)
            else:  # no crossover
                child1 = copy.deepcopy(parent1_transformed)
                child2 = copy.deepcopy(parent2_transformed)
                new_population += [child1, child2]

        new_population = new_population[:size]  # sbx_crossover creates pairs so this is for odd size of the population

        cauchy_mutation_population(new_population,
                                   bounds=self._bounds_transformed,
                                   gamma=self.__gamma_default,
                                   mutation_rate=self._mutation_rate,
                                   rng=self._rng)

        new_population = [self._transform_solution_back(sol) for sol in new_population]

        return new_population

    @staticmethod
    def get_elites(population, size=1) -> list:
        if not isinstance(size, int):
            raise TypeError
        if not (0 <= size <= len(population)):
            raise ValueError
        elites = sorted(population)[:size]
        return elites

    @staticmethod
    def update_population(population) -> None:
        for sol in population:
            sol.update()

    def is_solution_inside_bounds(self, sol, bounds=None) -> bool:
        if bounds is None:
            bounds = self._bounds
        return is_values_inside_bounds(sol.x, bounds)

    def filter_population(self, population) -> list:

        logger.info('filter_population: START')

        population_filtered = []
        for i, sol in enumerate(population):
            if not sol.is_updated():
                logger.info(f'  {i}: not updated')
                continue
            if not sol.is_valid():
                logger.info(f'  {i} not valid')
                continue
            if not self.is_solution_inside_bounds(sol):
                logger.info(f'  {i} outside bounds')
                continue

            logger.info(f'  {i} kept')
            population_filtered.append(sol)

        logger.info('filter_population: END')

        return population_filtered

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

        population = self.generate_population(n_solutions)
        for i in range(n_epochs):
            self.update_population(population)
            population = self.filter_population(population)
            elites = self.get_elites(population, n_elites)
            mutants = self.get_mutants(population, n_solutions - n_elites)
            population = elites + mutants
        return population
