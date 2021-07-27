import numpy as np
import copy

from ...helpers import uniform_vector


def cauchy_inverse_cdf(gamma, rng):
    if rng is None:
        rng = np.random.default_rng()
    return gamma * np.tan(np.pi * (rng.random() - 0.5))


def cauchy_mutation(genes, gamma=1, bounds=None, rng=None):  # do not change gamma=1

    if bounds is None:
        bounds = [[None, None]] * len(genes)
    assert (len(genes) == len(bounds))

    genes_new = []
    # if rng is None:
    #     rng = np.random.default_rng()
    shift = cauchy_inverse_cdf(gamma, rng)
    shift_vec = shift * uniform_vector(len(genes), rng=rng)  # vector mutation

    for gene, (lb, ub), shift in zip(genes, bounds, shift_vec):
        if lb is not None:  # bounce, TODO: rewrite for general case
            assert (lb <= gene <= ub), f"Violated: {lb} <= {gene} <= {ub}"
            b = ub - gene
            ptp = ub - lb
            while shift < 0:
                shift += 2 * ptp
            shift %= (2 * ptp)
            shift = np.abs(np.abs(shift - b) - ptp) - (ptp - b)

        gene_new = gene + shift
        genes_new.append(gene_new)

    return genes_new


def cauchy_mutation_population(population, bounds, gamma, mutation_rate, rng=None):

    if rng is None:
        rng = np.random.default_rng()

    if len(population):

        n_genes = len(population[0])
        if n_genes != len(bounds):
            raise ValueError

        genes = np.concatenate([organism.x.flatten() for organism in population])

        p = rng.random(len(population))
        shifts = gamma * np.tan(np.pi * (p - 0.5))

        mut_mask = rng.random(len(population)) < mutation_rate
        shifts = shifts * mut_mask

        shifts = np.tile(shifts, (n_genes, 1)).T.flatten()

        u = rng.standard_normal(n_genes * len(population)).reshape((n_genes, len(population)))
        u = u / np.linalg.norm(u, axis=1)[:, None]
        u = u.flatten()

        shifts = shifts * u
        lb, ub = np.tile(bounds, (len(population), 1)).T

        assert (np.all(lb <= genes) and np.all(genes <= ub)), "genes are outside bounds"

        ptp = ub - lb
        b = ub - genes

        shifts = np.remainder(shifts, 2 * ptp)
        shifts = np.abs(np.abs(shifts - b) - ptp) - (ptp - b)

        genes = genes + shifts
        genes = np.reshape(genes, (len(population), n_genes))

    else:
        genes = []

    mutants = copy.deepcopy(population)

    for i in range(len(mutants)):
        mutants[i].x = genes[i]

    return mutants


def reflection(ub, lb, genes, shifts):
    ptp = ub - lb
    b = ub - genes
    shifts = np.remainder(shifts, 2 * ptp)
    shifts = np.abs(np.abs(shifts - b) - ptp) - (ptp - b)
    genes = genes + shifts
    return genes
