#!/usr/bin/env python

import logging

from mpi4py import MPI
from tqdm.auto import tqdm
import numpy as np
from pypoptim.algorythm.ga import GA
from pypoptim.helpers import argmin, is_values_inside_bounds
from mpi_utils import allocate_recvbuf, allgather, population_from_recvbuf
from pypoptim.algorythm import Solution
from scipy.optimize import rosen


def error(x):
    return rosen(x)


class Sol(Solution):

    def update(self):
        self._y = error(self.x)

    def is_valid(self):
        return self.is_updated()


def mpi_script():
    logger = logging.getLogger(__name__)

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()

    n_generations = 1000
    n_genes = 20
    n_organisms_per_process = 100
    n_organisms = n_organisms_per_process * comm_size
    n_elites = 1 * comm_size

    bounds = np.tile([-5, 5], (n_genes, 1))

    recvbuf_dict = allocate_recvbuf(n_genes, n_organisms_per_process, comm)

    rng = np.random.Generator(np.random.PCG64(42 + comm_rank))
    ga_optim = GA(Sol,
                  bounds=bounds,
                  mutation_rate=0.1,
                  crossover_rate=0.9,
                  selection_force=2,
                  rng=rng)

    batch = ga_optim.generate_population(n_organisms_per_process)

    if comm_rank == 0:
        print(f'# comm_size = {comm_size}')
        pbar = tqdm(total=n_generations, ascii=True)

    loss_global = np.inf

    for epoch in range(n_generations):

        if comm_rank == 0:
            pbar.set_postfix_str("CALC")
        for i, sol in enumerate(batch):
            sol.update()

        if comm_rank == 0:
            pbar.set_postfix_str("GATHER")

        allgather(batch, recvbuf_dict, comm)
        population = population_from_recvbuf(recvbuf_dict, Sol, n_organisms, n_genes)

        shift = comm_rank * n_organisms_per_process
        assert all(np.all(sol_b.x == sol_p.x) for sol_b, sol_p in zip(batch, population[shift:]))

        if comm_rank == 0:
            pbar.set_postfix_str("SAVE")

        index_best = argmin(population)
        assert population[index_best] is min(population)
        comm_rank_best = index_best // n_organisms_per_process
        index_best_batch = index_best % n_organisms_per_process

        loss_current = min(population).y
        assert loss_current <= loss_global
        loss_global = loss_current
        if comm_rank == 0:
            pbar.set_description(f'Loss: {loss_global:.6f}')

        if comm_rank == comm_rank_best:
            sol_best = batch[index_best_batch]

            assert sol_best is min(batch)
            assert np.all(sol_best.x == min(population).x)
            assert sol_best.is_updated()
            assert sol_best.is_valid()
            assert is_values_inside_bounds(sol_best.x, bounds)
            assert ga_optim.is_solution_inside_bounds(sol_best)

        if comm_rank == 0:
            pbar.set_postfix_str("GENE")

        population = ga_optim.filter_population(population)

        population.sort()

        if len(population) <= 3:
            if comm_rank == 0:
                msg = f"# Not enough organisms for genetic operations left: {len(population)}"
                raise RuntimeError(msg)

        elites_all = population[:n_elites]  # len may be less than config['n_elites'] due to invalids
        elites_batch = elites_all[comm_rank::comm_size]  # elites_batch may be empty
        n_elites = len(elites_batch)
        n_mutants = n_organisms_per_process - n_elites

        mutants_batch = ga_optim.get_mutants(population, n_mutants)
        batch = elites_batch + mutants_batch

        assert (len(batch) == n_organisms_per_process)

        if comm_rank == 0:
            pbar.update(1)
            pbar.refresh()

    if comm_rank == 0:
        pbar.set_postfix_str("DONE")
        pbar.refresh()


if __name__ == '__main__':
    level = logging.INFO
    logging.basicConfig(level=level)
    logging.getLogger('numba').setLevel(logging.CRITICAL)  # https://stackoverflow.com/a/63471108/13213091
    mpi_script()
