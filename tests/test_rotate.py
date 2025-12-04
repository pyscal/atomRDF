import pytest
import os
from atomrdf import KnowledgeGraph, System
import numpy as np
import atomrdf.build as build
from atomrdf.transform import rotate

def test_rotate():
    kg = KnowledgeGraph()
    s = build.bulk("Al", graph=kg)
    rotation_vectors = [[ 1, 1, 0],
                [-1, 1, 0],
                [ 0, 0, 1]]
    srot = rotate(s, rotation_vectors)
    assert np.abs(srot.box[0][0] - 5.728) < 1E-3
    assert np.abs(srot.box[1][1] - 5.728) < 1E-3
    assert np.abs(srot.box[2][2] - 4.050) < 1E-3
    assert srot.natoms == 8

