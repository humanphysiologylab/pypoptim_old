import numpy as np
from sklearn.metrics import mean_squared_error as MSE


def RMSE(x, y):
    return MSE(x, y, squared=False)

def calculate_RMSE(x, y):
    assert(len(x) == len(y))
    return np.sqrt(np.mean(((x - y)**2)))

def calculate_RMSE_balanced(x, y):
    assert(len(x) == len(y))
    x = (x - y.min(axis=0)) / y.ptp(axis=0)
    y = (y - y.min(axis=0)) / y.ptp(axis=0)
    return np.sqrt(np.sum(((x - y)**2)) / len(x))

def calculate_composite_RMSE_V_CaT(x, y):
    # x -- model, y -- experiment
    assert (len(x) == len(y))

    v_model, cat_model = x.T
    v_exp, cat_exp = y.T

    cat_model = (cat_model - cat_model.min(axis=0)) / cat_model.ptp(axis=0)  # to balance V and CaT

    rmse_v = calculate_RMSE_balanced(v_model, v_exp)    # v_exp --> [0, 1]
    cat_exp_scaled, rmse_cat, coeffs = autoscaling(signal_to_scale=cat_exp,
                                                   signal_reference=cat_model)

    rmse_total = rmse_v + rmse_cat

    return rmse_total

def calculate_RMSE_weightened(x, y, weights):
    return np.sum(MSE(x, y, squared=False, multioutput='raw_values') * weights)


def autoscaling(signal_to_scale, signal_reference):
    def scalar_multiplications(a, b):  # TODO: REWRITE!!!
        assert (len(a) == len(b))
        assert (a.size == b.size)
        coefficients = np.array([np.dot(a, b), np.sum(a), np.sum(b), np.sum(a ** 2), len(a)])
        return coefficients

    c = scalar_multiplications(signal_to_scale, signal_reference)

    if c[1] == 0 or c[1] * c[1] - c[4] * c[3] == 0:
        alpha = 0
        beta = 0
    else:
        beta = (c[0] * c[1] - c[2] * c[3]) / (c[1] * c[1] - c[4] * c[3])
        alpha = (c[2] - beta * c[4]) / c[1]

    signal_scaled = signal_to_scale * alpha + beta
    rmse = RMSE(signal_scaled, signal_reference)  # TODO: RMSE should be removed from this function

    return signal_scaled, rmse, (alpha, beta)