import pytest
import os
from atomrdf import KnowledgeGraph
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.namespace import CMSO, PLDO, ASMO, CDCO
import shutil

import atomrdf.build as build


def test_structuregraph():
    s = KnowledgeGraph()
    sys = build.bulk("Fe", graph=s)

    vis = s.visualise()
    assert vis != None

    vis = s.visualise()
    assert vis != None

    s.write("temp.ttl", format="turtle")
    s = KnowledgeGraph(graph_file="temp.ttl")

    sys = build.defect.vacancy("Fe", cubic=True, graph=s, no_of_vacancies=1)

    s = KnowledgeGraph()
    sys = build.bulk("Fe", graph=s)
    assert s.n_samples == 1
    # res = s.query_sample("NumberOfAtoms", 2)
    # assert(len(res) == 1)


def test_logger():
    if os.path.exists("tests/atomrdf.log"):
        os.remove("tests/atomrdf.log")
    s = KnowledgeGraph(enable_log=True)
    s.log("testing-logger")
    assert str(type(s.log).__name__) == "method"


def test_add_cross_triple():
    s = KnowledgeGraph(enable_log=True)
    atoms = build.bulk("Fe", graph=s)
    sample_id = atoms.info["id"]
    # Get material URIRef from graph
    material = s.value(sample_id, CMSO.hasMaterial)
    status, _ = s._check_domain_if_uriref(
        (material, CDCO.hasCrystallographicDefect, PLDO.AntiphaseBoundary)
    )
    assert status == True


def test_add_quantity():
    """Test adding calculated quantities using pydantic model fields."""
    s = KnowledgeGraph(enable_log=True)
    atoms = build.bulk("Fe", graph=s)
    sample_id = atoms.info["id"]

    # Retrieve sample from graph
    sample = s.get_sample_as_structure(sample_id)
    assert sample is not None
    assert sample.material.crystal_structure.spacegroup_symbol == "Im-3m"


def test_archive():
    s = KnowledgeGraph(enable_log=True)
    sys = build.bulk("Fe", graph=s)
    sys = build.bulk("Cu", graph=s)
    if os.path.exists("test_archive.tar.gz"):
        os.remove("test_archive.tar.gz")
    if os.path.exists("test_archive"):
        shutil.rmtree("test_archive")
    s.archive("test_archive")
    assert os.path.exists("test_archive.tar.gz")

    s = KnowledgeGraph.unarchive("test_archive.tar.gz")
    assert s.n_samples == 2
    os.remove("test_archive.tar.gz")
    shutil.rmtree("test_archive")


def test_sparql_query():
    """Test SPARQL queries using graph.query_rdflib()."""
    kg = KnowledgeGraph()
    struct_Fe = build.bulk("Fe", cubic=True, graph=kg)  # bcc cubic cell, 2 atoms
    struct_Si = build.bulk("Si", cubic=True, graph=kg)  # cubic diamond, 8 atoms
    struct_Al = build.bulk("Al", cubic=True, graph=kg)  # fcc cubic cell, 4 atoms
    query = """
	PREFIX cmso: <http://purls.helmholtz-metadaten.de/cmso/>
	SELECT DISTINCT ?symbol
	WHERE {
	    ?sample cmso:hasNumberOfAtoms ?number .
	    ?sample cmso:hasMaterial ?material .
	    ?material cmso:hasStructure ?structure .
	    ?structure cmso:hasSpaceGroupSymbol ?symbol .
	FILTER (?number="4"^^xsd:integer)
	}"""
    # Use rdflib query directly via kg.graph.query()
    res = kg.graph.query(query)
    results = list(res)
    assert len(results) > 0
    # Al fcc cubic cell has 4 atoms and Fm-3m space group
    assert results[0][0].toPython() == "Fm-3m"


def test_extract_sample():
    """Test extracting sample as AtomicScaleSample structure."""
    kg = KnowledgeGraph()
    atoms = build.bulk("Fe", cubic=True, graph=kg)  # cubic cell gives 2 atoms for bcc
    sample_id = atoms.info["id"]

    # Use get_sample_as_structure to retrieve AtomicScaleSample
    sample = kg.get_sample_as_structure(sample_id)
    assert sample is not None
    # Convert to ASE Atoms to check positions
    atoms_from_sample = sample.to_structure()
    assert len(atoms_from_sample.positions) == 2
    assert sample.material.element_ratio["Fe"] == 1.0

    # Test to_file with default format (lammps-data)
    kg.to_file(sample_id, filename="POSCAR")
    assert os.path.exists("POSCAR")
    os.remove("POSCAR")

    # Test to_file with custom format
    kg.to_file(sample_id, filename="POSCAR", format="cif")
    assert os.path.exists("POSCAR")
    os.remove("POSCAR")


# def test_add_domain_ontoterm():
# 	from atomrdf.namespace import CMSO, PLDO
# 	s = KnowledgeGraph()
# 	sys = System.create.element.Fe(graph=s)
# 	status, _ = s._check_domain_if_ontoterm((CMSO.Material, CMSO.hasDefect, PLDO.AntiphaseBoundary))
# 	assert status == True


def test_purge():
    s = KnowledgeGraph()
    sys = build.bulk("Fe", graph=s)
    s.purge(force=True)
    assert s.n_samples == 0

    s = KnowledgeGraph(store="db", store_file=f"testr.db")
    sys = build.bulk("Fe", graph=s)
    s.purge(force=True)
    assert s.n_samples == 0
