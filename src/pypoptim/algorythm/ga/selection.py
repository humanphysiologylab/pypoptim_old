import random

def tournament_selection(p, k=2):
    return min(random.sample(p, k))
