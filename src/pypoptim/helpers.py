import itertools
import numpy as np
from numba import njit


def uniform_vector(n=1, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    u = rng.standard_normal(n)
    return u / np.linalg.norm(u)


def argmax(l):
    return max(enumerate(l), key=lambda x: x[1])[0]


def argmin(l):
    return min(enumerate(l), key=lambda x: x[1])[0]


def argmax_list_of_dicts(l, key):
    return max(enumerate(l), key=lambda x: x[1][key])[0]


def find_index_first(seq, condition):
    #  https://stackoverflow.com/a/8534381/13213091
    return next((i for i, x in enumerate(seq) if condition(x)), None)


def get_value_by_key(array, legend, key):
    return dict(zip(legend['name'], array))[key]


def flatten_iterable(x):
    return list(itertools.chain(*x))


def batches_from_list(l, n_batches=1):
    return [l[i::n_batches] for i in range(n_batches)]


def random_value_from_bounds(bounds, log_scale=False, rng=None):
    if len(bounds) != 2 or bounds[0] >= bounds[1]:
        raise ValueError
    if rng is None:
        rng = np.random.default_rng()
    r = rng.random()
    if log_scale:
        return (bounds[0] ** (1 - r)) * (bounds[1] ** r)
    else:
        return bounds[0] * (1 - r) + bounds[1] * r


def rastrigin(x, A=10):
    x = np.array(x)
    return sum(x**2 + A * (1 - np.cos(2 * np.pi * x)))


def ma(x, n):
    return np.convolve(x, np.ones(n)/n, mode='valid')


def calculate_mean_abs_noise(array, n=31):
    array_ma = np.apply_along_axis(func1d=ma, axis=0, arr=array, n=n)
    array_valid = array[n//2: n//2 + len(array_ma)]
    noise = np.mean(np.abs(array_valid - array_ma), axis=0)
    return noise


@njit
def transform_genes_bounds(genes, bounds, gammas, mask_multipliers):
    assert (len(genes) == len(bounds) == len(gammas) == len(mask_multipliers))

    genes_transformed = np.zeros_like(genes)
    bounds_transformed = np.zeros_like(bounds)

    scaler_dimensional = 1 / np.sqrt(len(genes))
    for i in range(len(genes)):
        lb, ub = bounds[i]
        gene = genes[i]
        if mask_multipliers[i]:  # log10 scale
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


@njit
def transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds_back, mask_multipliers):
    assert (len(genes_transformed) == len(bounds_transformed) == len(mask_multipliers))

    genes_back = np.zeros_like(genes_transformed)

    for i in range(len(genes_transformed)):  # log10 scale
        lb_back, ub_back = bounds_back[i]
        lb_tran, ub_tran = bounds_transformed[i]
        gene = genes_transformed[i]
        if mask_multipliers[i]:
            genes_back[i] = np.log10(lb_back) + (gene - lb_tran) / (ub_tran - lb_tran) * (
                    np.log10(ub_back) - np.log10(lb_back))
            genes_back[i] = np.power(10, genes_back[i])
        else:  # linear scale
            genes_back[i] = lb_back + (gene - lb_tran) / (ub_tran - lb_tran) * (ub_back - lb_back)

    return genes_back


def wrapped_ndptr(*args, **kwargs):
    #  https://stackoverflow.com/a/37664693/13213091
    base = np.ctypeslib.ndpointer(*args, **kwargs)
    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)
    return type(base.__name__, (base,), {'from_param': classmethod(from_param)})

DoubleArrayType_1D = wrapped_ndptr(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
DoubleArrayType_2D = wrapped_ndptr(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')


def strip_comments(code, comment_char='#'):
    lines = []
    for i, line in enumerate(code.split('\n')):
        items = line.split(comment_char)
        if len(items) >= 2:
            line = items[0]
        lines.append(line)
    code = '\n'.join(lines)
    return code


def autoscaling(signal_to_scale, signal_reference):

    def scalar_multiplications(a, b):
        if len(a) != len(b):
            raise ValueError
        coefficients = np.array([np.dot(a, b),
                                 np.sum(a),
                                 np.sum(b),
                                 np.sum(a**2),
                                 len(a)])
        return coefficients

    c = scalar_multiplications(signal_to_scale, signal_reference)

    if c[1] == 0 or c[1] * c[1] - c[4] * c[3] == 0:
        alpha = 0
        beta = 0
    else:
        beta = (c[0] * c[1] - c[2] * c[3]) / (c[1] * c[1] - c[4] * c[3])
        alpha = (c[2] - beta * c[4]) / c[1]

    signal_scaled = signal_to_scale * alpha + beta

    return signal_scaled, (alpha, beta)
