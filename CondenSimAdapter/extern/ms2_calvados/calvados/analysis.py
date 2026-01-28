"""
Minimal analysis module for ms2_calvados
Only keep self_distances function (used by components.py)
"""

import numpy as np
from MDAnalysis.analysis import distances


def self_distances(pos, box=None):
    """ Self distance map for matrix of positions

    If box dimensions are provided, distances are
    calculated using minimum image convention

    Input: Matrix of positions and (optional) box dimensions
    Output: Self distance map
    """
    N = len(pos)
    dmap = np.zeros((N, N))
    if box is not None:
        d = distances.self_distance_array(pos, box)
    else:
        d = distances.self_distance_array(pos)
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            dmap[i, j] = d[k]
            dmap[j, i] = d[k]
            k += 1
    return dmap








