import numpy as np

from ..helpers import uniform_vector

def test_uniform_vector():

    n = 42
    uv = uniform_vector(n)
    
    assert len(uv) == n
    assert np.isclose(np.linalg.norm(uv), 1)
