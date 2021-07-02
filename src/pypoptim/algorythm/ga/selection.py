import random
from deprecated import deprecated


@deprecated
def tournament_selection(p, k=2, key='loss'):
    return min(random.sample(p, k), key=lambda x: x[key])
