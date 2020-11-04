from mpi4py import MPI
import os
import json
import copy
import time
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.integrate import odeint

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from src.model.maleckar import init_states_constants, compute_rates, legend
from src.helpers import get_value_by_key, update_array_from_kwargs,\
    calculate_RMSE, flatten_list, batches_from_list, value_from_bounds, argmax_list_of_dicts
from src.algorythm.ga import do_step


def generate_organism(config, legend):
    genes = []
    state = []

    for mult in config['multipliers']:
        value = value_from_bounds(mult['bounds'], log_scale=True)
        genes.append(value)
        name = mult['name']

    for exp_cond in config['experimental_conditions']:

        initial_state = np.copy(exp_cond['initial_state'])
        state.append(initial_state)

        if "params" in exp_cond:
            for param in exp_cond["params"]:
                if "bounds" in param:
                    if not ("value" in param):
                        if param['name'] in legend['states']['name'].values:
                            value = get_value_by_key(initial_state, legend['states'], param['name'])
                            genes.append(value)
                            name = param['name']
                        else:
                            raise Exception(f"Invalid config entry:\n\t{param}")

                    else:
                        raise Exception(f"Invalid config entry:\n\t{param}")

    organism = dict(genes=np.array(genes), state=np.array(state))
    return organism


def init_population(config, legend):
    population = [generate_organism(config, legend) for _ in range(config['n_organisms'])]
    return population


def run_model(S, C, R, A, config):

    stim_period = get_value_by_key(C, legend['constants'], 'stim_period')
    t_sampling = config['t_sampling']
    n_beats = config['n_beats']
    t_space = np.arange(0, n_beats * stim_period, t_sampling)
    kwargs_odeint = dict(hmax=t_sampling, full_output=False)
    res = odeint(func=compute_rates, y0=S, t=t_space, args=(C, R, A), tfirst=True, **kwargs_odeint)

    return res


def update_phenotype_state(organism, legend, config):

    organism['phenotype'] = []

    for i, exp_cond in enumerate(config['experimental_conditions']):

        multipliers_names = []
        params_names = []

        update_kwargs = OrderedDict()
        index_gene = 0
        for mult in config['multipliers']:
            name = mult['name']
            multipliers_names.append(name)
            value = organism['genes'][index_gene]
            update_kwargs[name] = value
            index_gene += 1

        if "params" in exp_cond:
            for param in exp_cond["params"]:
                if "bounds" in param:
                    if not ("value" in param):
                        if param['name'] in legend['states']['name'].values:
                            name = param['name']
                            params_names.append(name)
                            value = organism['genes'][index_gene]
                            update_kwargs[name] = value
                            index_gene += 1
                        else:
                            raise Exception(f"Invalid config entry:\n\t{param}")
                    else:
                        name = param['name']
                        value = param['value']
                        update_kwargs[name] = value
                        index_gene += 1

        S, C = init_states_constants()
        R = np.zeros_like(S)
        A = np.zeros(len(legend['algebraic']))

        C_default = pd.Series(data=C, copy=True,
                              index=legend['constants']['name'])
        for name in multipliers_names:
            update_kwargs[name] *= C_default[name]

        C = update_array_from_kwargs(C, legend['constants'], **update_kwargs)
        S = organism['state'][i]
        S = update_array_from_kwargs(S, legend['states'], **update_kwargs)

        CL = exp_cond['CL']
        stim_period = CL / 1000.
        C = update_array_from_kwargs(C, legend['constants'], stim_period=stim_period)

        res = run_model(S, C, R, A, config)

        t_sampling = config['t_sampling']
        organism['phenotype'].append(res[-int(stim_period / t_sampling):, 0])
        organism['state'][i] = res[-1]

        for j, name in enumerate(params_names):
            assert(name in legend['states']['name'].values)
            n_multipliers = len(multipliers_names)
            value = get_value_by_key(organism['state'][i], legend['states'], name)
            organism['genes'][n_multipliers + j] = value

    return #  copy.deepcopy(organism)


def update_fitness(organism, config):
    loss = 0
    for i, exp_cond in enumerate(config['experimental_conditions']):
        phenotype_control = exp_cond['phenotype']
        phenotype_organism = organism['phenotype'][i][:len(phenotype_control)]
        loss += calculate_RMSE(phenotype_control, phenotype_organism)
    organism['fitness'] = -loss
    return #  copy.deepcopy(organism)


def generate_bounds_gammas(config):
    bounds = []
    gammas = []

    for mult in config['multipliers']:
        bounds.append(mult['bounds'])
        gammas.append(mult['gamma'])

    for exp_cond in config['experimental_conditions']:
        if "params" in exp_cond:
            for param in exp_cond["params"]:
                condition = ("bounds" in param) and ("gamma" in param) and not ("value" in param)
                if condition:
                    bounds.append(param['bounds'])
                    gammas.append(param['gamma'])

    bounds, gammas = np.array(bounds), np.array(gammas)
    return bounds, gammas


def save_epoch(population, dump_filename, output_folder_name):
    population = sorted(population, key=lambda organism: organism['fitness'], reverse=True)
    genes = np.array([organism['genes'] for organism in population])
    fitness = np.array([organism['fitness'] for organism in population])

    dump_current = np.hstack([genes, fitness[:, None]])
    with open(dump_filename, 'ba+') as file_dump:
        dump_current.tofile(file_dump)

    print(np.array2string(dump_current[0, :], formatter={'float_kind': lambda x: "%.3f" % x}, max_line_width=1000)[1:-1], flush=True)


def remove_invalids(population, bounds):
    population_valid = []
    for organism in population:
        condition_invalid = np.any(np.isnan(organism['genes'])) \
                            or np.any(np.isnan(organism['state'])) \
                            or (organism['fitness'] == np.NINF) \
                            or np.any(organism['genes'] <= bounds[:, 0]) \
                            or np.any(organism['genes'] >= bounds[:, 1])
        if not condition_invalid:
            population_valid.append(organism)
    return population_valid


def transform_genes_bounds(genes, bounds, gammas, n_multipliers):
    assert (len(genes) == len(bounds) == len(gammas))

    genes_transformed = np.zeros_like(genes)
    bounds_transformed = np.zeros_like(bounds)

    scaler_dimensional = 1 / np.sqrt(len(genes))
    for i in range(len(genes)):
        lb, ub = bounds[i]
        gene = genes[i]
        if i < n_multipliers:  # log10 scale
            bounds_transformed[i, 1] = np.log10(ub / lb) * 1 / (gammas[i] / scaler_dimensional)
            genes_transformed[i] = np.log10(gene)
            lb_temp = np.log10(lb)
            ub_temp = np.log10(ub)
        else:  # linear scale
            genes_transformed[i] = gene
            bounds_transformed[i, 1] = (ub - lb) * 1 / (gammas[i] / scaler_dimensional)
            lb_temp = lb
            ub_temp = ub
        genes_transformed[i] = (genes_transformed[i] - lb_temp) / (ub_temp - lb_temp) * bounds_transformed[i, 1]

    return genes_transformed, bounds_transformed


def transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds_back, n_multipliers):
    assert (len(genes_transformed) == len(bounds_transformed))

    genes_back = np.zeros_like(genes_transformed)

    for i in range(n_multipliers):  #log10 scale
        lb_back, ub_back = bounds_back[i]
        lb_tran, ub_tran = bounds_transformed[i]
        gene = genes_transformed[i]
        genes_back[i] = np.log10(lb_back) + (gene - lb_tran) / (ub_tran - lb_tran) * (
                    np.log10(ub_back) - np.log10(lb_back))
        genes_back[i] = np.power(10, genes_back[i])
    for i in range(n_multipliers, len(genes_transformed)):  # linear scale
        lb_back, ub_back = bounds_back[i]
        lb_tran, ub_tran = bounds_transformed[i]
        gene = genes_transformed[i]
        genes_back[i] = lb_back + (gene - lb_tran) / (ub_tran - lb_tran) * (ub_back - lb_back)

    return genes_back

# Initializing

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

if len(sys.argv) < 2:
    if comm_rank == 0:
        print(f"Usage: mpiexec -n 2 python {sys.argv[0]} config.json")
    exit()

config_filename = sys.argv[1]
config_path = os.path.dirname(os.path.realpath(config_filename))

with open(config_filename) as f:
    config = json.load(f)

if config['n_organisms'] % comm_size != 0:
    n_organisms_new = int(np.ceil(config['n_organisms'] / comm_size) * comm_size)
    if comm_rank == 0:
        print(f'# `n_organisms` is changed from {config["n_organisms"]} to {n_organisms_new}', flush=True)
    config['n_organisms'] = n_organisms_new


for exp_cond in config['experimental_conditions']:
    CL = exp_cond['CL']
    filename_phenotype = os.path.join(config_path, exp_cond['filename_phenotype'])
    exp_cond['phenotype'] = np.loadtxt(filename_phenotype)
    filename_state = os.path.join(config_path, exp_cond['filename_state'])
    exp_cond['initial_state'] = np.loadtxt(filename_state)

output_folder_name = "./output"
output_folder_name_phenotype = os.path.join(output_folder_name, "phenotype")
output_folder_name_state = os.path.join(output_folder_name, "state")
dump_filename = os.path.join(output_folder_name, "dump")

if comm_rank == 0:
    for folder in output_folder_name_phenotype, output_folder_name_state:
        os.makedirs(folder, exist_ok=True)
    with open(dump_filename, "wb") as file_dump:  # create or clear and close
        pass
    time_start = time.time()
    print("# SIZE =", comm_size, flush=True)

bounds, gammas = generate_bounds_gammas(config)

n_orgsnisms_per_process = config['n_organisms'] // comm_size

organism_dummy = generate_organism(config, legend)
update_phenotype_state(organism_dummy, legend, config)

genes_size = organism_dummy['genes'].size
genes_shape = organism_dummy['genes'].shape

state_size = organism_dummy['state'].size
state_shape = organism_dummy['state'].shape

phenotype_size = np.concatenate(organism_dummy['phenotype']).size
phenotype_lens = [len(phenotype) for phenotype in organism_dummy['phenotype']]
phenotype_split_indices = np.cumsum(phenotype_lens)[:-1]

del organism_dummy

#  Main cycle
times = dict()

if comm_rank == 0:
    population = init_population(config=config, legend=legend)
    population = batches_from_list(population, comm_size)
else:
    population = None

batch = comm.scatter(population, root=0)

for epoch in range(config['n_generations']):

    #  calculations
    times['calc'] = time.time()
    for organism in batch:
        try:
            update_phenotype_state(organism, legend, config)
        except:
            print("Bad guy")
            organism['phenotype'] = [np.empty(length) for length in phenotype_lens]  # needed for gather
            organism['fitness'] = np.NINF
            continue
        update_fitness(organism, config)
    times['calc'] = time.time() - times['calc']

    times['gather_concat'] = time.time()
    sendbuf_genes = np.concatenate([organism['genes'].flatten() for organism in batch])
    sendbuf_state = np.concatenate([organism['state'].flatten() for organism in batch])
    #sendbuf_phenotype = np.concatenate([np.concatenate(organism['phenotype']) for organism in batch])
    sendbuf_fitness = np.array([organism['fitness'] for organism in batch])
    times['gather_concat'] = time.time() - times['gather_concat']

    times['gather_alloc'] = time.time()
    recvbuf_genes = np.empty([comm_size, n_orgsnisms_per_process * genes_size])
    recvbuf_state = np.empty([comm_size, n_orgsnisms_per_process * state_size])
    #recvbuf_phenotype = np.empty([comm_size, n_orgsnisms_per_process * phenotype_size])
    recvbuf_fitness = np.empty([comm_size, n_orgsnisms_per_process * 1])
    times['gather_alloc'] = time.time() - times['gather_alloc']

    times['gather'] = time.time()
    comm.Allgatherv(sendbuf_genes, recvbuf_genes)
    comm.Allgatherv(sendbuf_state, recvbuf_state)
    #comm.Allgatherv(sendbuf_phenotype, recvbuf_phenotype)
    comm.Allgatherv(sendbuf_fitness, recvbuf_fitness)
    times['gather'] = time.time() - times['gather']

    times['gather_reshape'] = time.time()
    recvbuf_genes = recvbuf_genes.reshape((config['n_organisms'], genes_size))
    recvbuf_state = recvbuf_state.reshape((config['n_organisms'], *state_shape))
    #recvbuf_phenotype = recvbuf_phenotype.reshape((config['n_organisms'], phenotype_size))
    recvbuf_fitness = recvbuf_fitness.flatten()
    times['gather_reshape'] = time.time() - times['gather_reshape']

    times['gather_pop'] = time.time()
    population = [dict(genes     = recvbuf_genes[i].copy(),
                       state     = recvbuf_state[i].copy(),
                       #phenotype = np.split(recvbuf_phenotype[i].copy(), indices_or_sections=phenotype_split_indices),
                       fitness   = recvbuf_fitness[i].copy()) for i in range(config['n_organisms'])]
    times['gather_pop'] = time.time() - times['gather_pop']

    index_best = argmax_list_of_dicts(population, 'fitness')
    comm_rank_best = index_best // n_orgsnisms_per_process
    index_best_batch = index_best % n_orgsnisms_per_process

    if comm_rank == comm_rank_best:
        organism_best = batch[index_best_batch]
        for i, exp_cond in enumerate(config['experimental_conditions']):
            CL = exp_cond['CL']
            np.savetxt(os.path.join(output_folder_name_phenotype, f"phenotype{CL}_{epoch}.txt"), organism_best['phenotype'][i])
            np.savetxt(os.path.join(output_folder_name_state, f"state_{CL}_{epoch}.txt"), organism_best['state'][i])

    if comm_rank == 0:
        save_epoch(population, dump_filename, output_folder_name)

    #  Next generation
    times['genetic'] = time.time()
    population = remove_invalids(population, bounds)

    assert(len(population) > 3), "Not enough organisms for genetic operations"

    population.sort(key=lambda organism: organism['fitness'], reverse=True)
    elites_all = population[:config['n_elites']]
    elites_batch = elites_all[comm_rank * n_orgsnisms_per_process: (comm_rank + 1) * n_orgsnisms_per_process]
    n_elites = len(elites_batch)
    #  elites_batch may be empty list

    for organism in population:
        organism['genes'], bounds_transformed = transform_genes_bounds(organism['genes'], bounds, gammas,
                                                                       n_multipliers=len(config['multipliers']))

    batch = do_step(population, new_size=n_orgsnisms_per_process - n_elites,
                    elite_size=0, bounds=bounds_transformed)

    batch += elites_batch

    for organism in batch:
        organism['genes'] = transform_genes_bounds_back(organism['genes'], bounds_transformed, bounds_back=bounds,
                                                        n_multipliers=len(config['multipliers']))
    times['genetic'] = time.time() - times['genetic']

    times['total'] = sum(times[name] for name in times if name != 'total')

    if comm_rank == 0:
        print(f"# EPOCH {epoch}:")
        for t in times:
            print(f"# {t}:\t{times[t]:.6f}\t{100 * times[t] / times['total']:.2f} %")

if comm_rank == 0:
    time_end = time.time()
    print("# TIME =", time_end - time_start, flush=True)
    file_dump.close()
