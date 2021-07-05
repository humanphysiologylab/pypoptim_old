import numpy as np

from ..helpers import uniform_vector
from ..helpers import transform_genes_bounds, transform_genes_bounds_back

def test_uniform_vector():

    n = 42
    uv = uniform_vector(n)

    assert len(uv) == n
    assert np.isclose(np.linalg.norm(uv), 1)


def test_transform_genes():

    genes  = np.array( [1.,      1.0,    -1.,   100.0])

    bounds = np.array([[0.,      0.1,    -3.,    10.0],
                       [2.,     10.0,     0.,  1000.0]]).T

    gammas = np.array( [1.,      1.0,     2.,     0.5])

    mask_multipliers = np.array([False, True, False, True])

    genes_transformed, bounds_transformed = transform_genes_bounds(genes, bounds, gammas, mask_multipliers)

    assert np.allclose(genes_transformed, np.array([0.5, 0.5, 0.5, 1]))

    assert np.allclose(bounds_transformed, np.array([[0.  , 1.  ],
                                                     [0.  , 1.  ],
                                                     [0.  , 0.75],
                                                     [0.  , 2.  ]]))

    genes_back = transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds, mask_multipliers)

    assert np.allclose(genes, genes_back)