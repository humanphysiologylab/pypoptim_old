import numpy as np


def tournament_selection(p, k=2, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    return min(rng.choice(p, k, replace=False))
