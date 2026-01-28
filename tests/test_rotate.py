import pytest
import os
from atomrdf import KnowledgeGraph
import numpy as np
import atomrdf.build as build
from atomrdf.transform import rotate


def test_rotate():
    """Test rotation of atomic structures."""
    kg = KnowledgeGraph()
    atoms = build.bulk("Al", cubic=True, graph=kg)

    # Define rotation vectors
    rotation_vectors = [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]

    # Rotate the structure
    atoms_rot = rotate(atoms, rotation_vectors)

    # Basic checks
    assert atoms_rot is not None
    assert len(atoms_rot) > 0

    # Check that cell exists
    cell = atoms_rot.get_cell()
    assert cell is not None
    assert len(cell) == 3

    # Verify the structure has been rotated (cell should be different)
    original_cell = atoms.get_cell()
    assert not np.allclose(cell, original_cell)
