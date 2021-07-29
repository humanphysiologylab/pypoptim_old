import numpy as np
import pandas as pd


def allocate_recvbuf(n_genes, n_organisms_per_process, comm):
    comm_size = comm.Get_size()

    # n_orgsnisms_per_process = config['runtime']['n_organisms'] // comm_size
    # config['runtime']['n_orgsnisms_per_process'] = n_orgsnisms_per_process

    genes_size = n_genes
    # config['runtime']['genes_size'] = genes_size

    # state_size = config['runtime']['states_initial'].size
    # state_shape = config['runtime']['states_initial'].shape
    # config['runtime']['state_shape'] = state_shape

    recvbuf_genes = np.empty([comm_size, n_organisms_per_process * genes_size])
    recvbuf_loss = np.empty([comm_size, n_organisms_per_process * 1])

    recvbuf_dict = dict(genes=recvbuf_genes,
                        loss=recvbuf_loss,
                        )

    return recvbuf_dict


def allgather(batch, recvbuf_dict, comm):
    sendbuf_genes = np.concatenate([sol.x for sol in batch])
    # sendbuf_state = np.concatenate([sol['state'].values.flatten() for sol in batch])
    sendbuf_loss = np.array([sol.y for sol in batch])
    # sendbuf_status = np.array([sol.status for sol in batch]).astype(float)

    comm.Allgatherv(sendbuf_genes, recvbuf_dict['genes'])
    # comm.Allgatherv(sendbuf_state,  recvbuf_dict['state'])
    comm.Allgatherv(sendbuf_loss, recvbuf_dict['loss'])
    # comm.Allgatherv(sendbuf_status, recvbuf_dict['status'])


def population_from_recvbuf(recvbuf_dict, SolModel, n_organisms, n_genes):
    recvbuf_genes = recvbuf_dict['genes']
    recvbuf_loss = recvbuf_dict['loss']

    recvbuf_genes = recvbuf_genes.reshape((n_organisms, n_genes))
    recvbuf_loss = recvbuf_loss.flatten()

    population = []

    for i in range(n_organisms):
        sol = SolModel(recvbuf_genes[i].copy())
        sol._y = recvbuf_loss[i].copy()
        population.append(sol)

    return population
