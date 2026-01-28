import pytest
from atomrdf import KnowledgeGraph
from atomrdf.datamodels.structure import AtomicScaleSample
import atomrdf.build as build


def test_structuregraph():
    """Test basic structure graph creation and operations."""
    s = KnowledgeGraph()

    # Test 1: Create simple bulk structure
    atoms = build.bulk("Fe", graph=s)
    assert atoms.info["id"] is not None
    assert s.n_samples == 1

    # Test 2: Create another structure with explicit parameters
    atoms = build.bulk("Cu", cubic=True, graph=s)
    assert atoms.info["id"] is not None
    assert s.n_samples == 2

    # Test 3: Read from file if it exists
    try:
        sample = AtomicScaleSample.from_file(
            "tests/al_data/Al.poscar", format="vasp", graph=s
        )
        assert sample.id is not None
        assert s.n_samples == 3
    except FileNotFoundError:
        # File doesn't exist, skip this part
        pass
