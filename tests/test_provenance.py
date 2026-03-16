"""Tests for atomrdf.io.provenance – Provenance tracing module."""

import pytest
from ase import Atoms

from atomrdf import KnowledgeGraph
from atomrdf.io import from_workflow_input, Provenance
from atomrdf.namespace import ASMO, PROV


# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #


def _two_step_workflow_data():
    """Return a dict with 3 samples and 2 workflow steps chained together.

    Chain: sample_A → sim1 → sample_B → sim2 → sample_C
    sim2 has a calculated property.
    """
    return {
        "computational_sample": [
            {
                "id": "sample_A",
                "material": {
                    "element_ratio": {"Fe": 1.0},
                    "crystal_structure": {
                        "spacegroup_symbol": "Im-3m",
                        "spacegroup_number": 229,
                        "unit_cell": {
                            "bravais_lattice": "bcc",
                            "lattice_parameter": [2.87, 2.87, 2.87],
                            "angle": [90, 90, 90],
                        },
                    },
                },
                "simulation_cell": {
                    "volume": {"value": 23.55},
                    "number_of_atoms": 2,
                    "length": [2.87, 2.87, 2.87],
                    "vector": [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
                    "angle": [90, 90, 90],
                },
                "atom_attribute": {
                    "position": [[0, 0, 0], [1.435, 1.435, 1.435]],
                    "species": ["Fe", "Fe"],
                },
            },
            {
                "id": "sample_B",
                "material": {
                    "element_ratio": {"Fe": 1.0},
                    "crystal_structure": {
                        "spacegroup_symbol": "Im-3m",
                        "spacegroup_number": 229,
                        "unit_cell": {
                            "bravais_lattice": "bcc",
                            "lattice_parameter": [2.87, 2.87, 2.87],
                            "angle": [90, 90, 90],
                        },
                    },
                },
                "simulation_cell": {
                    "volume": {"value": 23.55},
                    "number_of_atoms": 2,
                    "length": [2.87, 2.87, 2.87],
                    "vector": [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
                    "angle": [90, 90, 90],
                },
                "atom_attribute": {
                    "position": [[0.01, 0.01, 0.01], [1.44, 1.44, 1.44]],
                    "species": ["Fe", "Fe"],
                },
            },
            {
                "id": "sample_C",
                "material": {
                    "element_ratio": {"Fe": 1.0},
                    "crystal_structure": {
                        "spacegroup_symbol": "Im-3m",
                        "spacegroup_number": 229,
                        "unit_cell": {
                            "bravais_lattice": "bcc",
                            "lattice_parameter": [2.87, 2.87, 2.87],
                            "angle": [90, 90, 90],
                        },
                    },
                },
                "simulation_cell": {
                    "volume": {"value": 23.55},
                    "number_of_atoms": 2,
                    "length": [2.87, 2.87, 2.87],
                    "vector": [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
                    "angle": [90, 90, 90],
                },
                "atom_attribute": {
                    "position": [[0.02, 0.02, 0.02], [1.45, 1.45, 1.45]],
                    "species": ["Fe", "Fe"],
                },
            },
        ],
        "workflow": [
            {
                "method": "MolecularStatics",
                "input_sample": ["sample_A"],
                "output_sample": ["sample_B"],
                "degrees_of_freedom": ["AtomicPositionRelaxation"],
                "interatomic_potential": {
                    "potential_type": "EmbeddedAtomModel",
                    "uri": "https://doi.org/10.xxxx",
                },
                "software": {
                    "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
                    "label": "LAMMPS",
                },
            },
            {
                "method": "MolecularStatics",
                "algorithm": "EquationOfStateFit",
                "input_sample": ["sample_B"],
                "output_sample": ["sample_C"],
                "calculated_property": [
                    {
                        "label": "EquilibriumEnergy",
                        "value": -4.316,
                        "unit": "EV",
                        "associate_to_sample": ["sample_C"],
                        "basename": "TotalEnergy",
                    },
                    {
                        "label": "BulkModulus",
                        "value": 1.3,
                        "unit": "GigaPA",
                        "associate_to_sample": ["sample_C"],
                        "basename": "BulkModulus",
                    },
                ],
                "degrees_of_freedom": [
                    "AtomicPositionRelaxation",
                    "CellVolumeRelaxation",
                ],
                "interatomic_potential": {
                    "potential_type": "EmbeddedAtomModel",
                    "uri": "https://doi.org/10.xxxx",
                },
                "software": {
                    "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
                    "label": "LAMMPS",
                },
            },
        ],
    }


@pytest.fixture
def two_step_kg():
    """Build a KG with a 2-step provenance chain and return (kg, result)."""
    kg = KnowledgeGraph()
    result = from_workflow_input(_two_step_workflow_data(), graph=kg)
    return kg, result


# ------------------------------------------------------------------ #
# Tests – Provenance.from_sample                                       #
# ------------------------------------------------------------------ #


def test_from_sample_basic(two_step_kg):
    """Tracing from the final sample should recover 2 steps."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    assert len(prov) == 2
    assert prov[0]["activity_type"] in ("Simulation", "EnergyCalculation")
    assert prov[1]["activity_type"] in ("Simulation", "EnergyCalculation")


def test_from_sample_order(two_step_kg):
    """Steps should be root-first: sim1 before sim2."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]
    sample_a_uri = result["sample_map"]["sample_A"]
    sample_b_uri = result["sample_map"]["sample_B"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    # First step: A → B
    assert prov[0]["input_sample_id"] == sample_a_uri
    assert prov[0]["output_sample_id"] == sample_b_uri

    # Second step: B → C
    assert prov[1]["input_sample_id"] == sample_b_uri
    assert prov[1]["output_sample_id"] == sample_c_uri


def test_from_sample_ase_structures(two_step_kg):
    """input_sample and output_sample should be ASE Atoms objects."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    for step in prov:
        assert isinstance(step["output_sample"], Atoms)
        if step["input_sample_id"] is not None:
            assert isinstance(step["input_sample"], Atoms)


def test_from_sample_single_step(two_step_kg):
    """Tracing from sample_B should find exactly 1 step."""
    kg, result = two_step_kg
    sample_b_uri = result["sample_map"]["sample_B"]

    prov = Provenance.from_sample(kg, sample_b_uri)
    assert len(prov) == 1


def test_from_sample_root(two_step_kg):
    """Tracing from sample_A (root) should find 0 steps."""
    kg, result = two_step_kg
    sample_a_uri = result["sample_map"]["sample_A"]

    prov = Provenance.from_sample(kg, sample_a_uri)
    assert len(prov) == 0


# ------------------------------------------------------------------ #
# Tests – Provenance.from_property                                     #
# ------------------------------------------------------------------ #


def test_from_property_basic(two_step_kg):
    """Tracing from a calculated property should also recover 2 steps."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    # Find a calculated property URI on sample_C
    from rdflib import URIRef

    props = [
        str(o)
        for _, _, o in kg.triples(
            (URIRef(sample_c_uri), ASMO.hasCalculatedProperty, None)
        )
    ]
    assert len(props) >= 1, "sample_C should have at least one calculated property"

    prov = Provenance.from_property(kg, props[0])
    assert len(prov) == 2


def test_from_property_raises_on_invalid(two_step_kg):
    """from_property should raise ValueError for a bogus URI."""
    kg, _ = two_step_kg
    with pytest.raises(ValueError, match="No sample found"):
        Provenance.from_property(kg, "http://bogus/property/xyz")


# ------------------------------------------------------------------ #
# Tests – step dict contents (simulation metadata)                     #
# ------------------------------------------------------------------ #


def test_step_method_populated(two_step_kg):
    """The method field should be populated from Simulation.from_graph."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    for step in prov:
        assert step["method"] is not None, "method should be populated"


def test_step_algorithm_on_second_step(two_step_kg):
    """The second step has EquationOfStateFit algorithm."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    # Second step should have algorithm
    assert prov[1]["algorithm"] is not None


def test_step_calculated_properties(two_step_kg):
    """The second step should have calculated properties."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    calc_props = prov[1]["calculated_properties"]
    assert len(calc_props) >= 1


def test_step_degrees_of_freedom(two_step_kg):
    """Steps should have degrees_of_freedom populated."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    assert len(prov[0]["degrees_of_freedom"]) >= 1
    assert len(prov[1]["degrees_of_freedom"]) >= 1


def test_step_software(two_step_kg):
    """Steps should have software populated."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    assert len(prov[0]["software"]) >= 1


def test_step_interatomic_potential(two_step_kg):
    """Steps should have interatomic_potential populated."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)

    assert prov[0]["interatomic_potential"] is not None


# ------------------------------------------------------------------ #
# Tests – kg.trace convenience method                                  #
# ------------------------------------------------------------------ #


def test_kg_trace_sample(two_step_kg):
    """kg.trace(sample_uri) should return a Provenance."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = kg.trace(sample_c_uri)
    assert isinstance(prov, Provenance)
    assert len(prov) == 2


def test_kg_trace_property(two_step_kg):
    """kg.trace(property_uri) should detect it's a property and trace."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    from rdflib import URIRef

    props = [
        str(o)
        for _, _, o in kg.triples(
            (URIRef(sample_c_uri), ASMO.hasCalculatedProperty, None)
        )
    ]
    assert len(props) >= 1

    prov = kg.trace(props[0])
    assert isinstance(prov, Provenance)
    assert len(prov) == 2


# ------------------------------------------------------------------ #
# Tests – visualize                                                    #
# ------------------------------------------------------------------ #


def test_visualize_returns_digraph(two_step_kg):
    """visualize() should return a graphviz.Digraph."""
    import graphviz

    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)
    dot = prov.visualize()

    assert isinstance(dot, graphviz.Digraph)


# ------------------------------------------------------------------ #
# Tests – repr / iteration                                             #
# ------------------------------------------------------------------ #


def test_repr(two_step_kg):
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)
    r = repr(prov)
    assert "2 steps" in r


def test_iterable(two_step_kg):
    """Provenance should be iterable and subscriptable."""
    kg, result = two_step_kg
    sample_c_uri = result["sample_map"]["sample_C"]

    prov = Provenance.from_sample(kg, sample_c_uri)
    steps = list(prov)
    assert len(steps) == 2
    assert prov[0] is steps[0]
