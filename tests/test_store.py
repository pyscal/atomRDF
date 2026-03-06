import pytest
import os
import shutil
import tempfile

import oxrdflib  # noqa: F401 — must be installed; CI installs via environment.yml

from atomrdf import KnowledgeGraph
from atomrdf.namespace import CMSO, PLDO
import atomrdf.build as build


def test_sqlalchemy():
    s = KnowledgeGraph(store="SQLAlchemy", store_file="aa")
    sys = build.bulk("Fe", graph=s)


# ─── Oxigraph store tests ────────────────────────────────────────────────────


@pytest.fixture()
def tmp_store_dir():
    """Temporary directory cleaned up after each test."""
    d = tempfile.mkdtemp(prefix="atomrdf_oxigraph_")
    yield d
    shutil.rmtree(d, ignore_errors=True)


def _build_sample_kg(store_kwargs: dict) -> KnowledgeGraph:
    """Helper: create a KG with the given store kwargs and add a bulk Fe cell."""
    kg = KnowledgeGraph(**store_kwargs)
    build.bulk("Fe", graph=kg)
    return kg


def test_oxigraph_inmemory_basic():
    """In-memory Oxigraph KG: can add triples and retrieve them."""
    kg = _build_sample_kg({"store": "Oxigraph"})
    n = len(kg.graph)
    assert n > 0, f"Expected triples in graph, got {n}"


def test_oxigraph_inmemory_sample_count():
    """In-memory Oxigraph KG: AtomicScaleSample is created."""
    kg = _build_sample_kg({"store": "Oxigraph"})
    assert kg.n_samples == 1


def test_oxigraph_inmemory_sparql():
    """In-memory Oxigraph KG: SPARQL SELECT returns results."""
    kg = _build_sample_kg({"store": "Oxigraph"})
    # Use the full purls IRI directly — CMSO.AtomicScaleSample is an OntoTerm,
    # not a URIRef, so str() gives the prefixed form.
    iri = "http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample"
    q = "SELECT (COUNT(?s) AS ?n) WHERE { ?s a <%s> . }" % iri
    rows = list(kg.graph.query(q))
    assert int(rows[0][0]) >= 1


def test_oxigraph_file_persistence(tmp_store_dir):
    """File-backed Oxigraph KG: triples survive close and reopen."""
    store_path = os.path.join(tmp_store_dir, "kg_persist")

    kg = _build_sample_kg({"store": "Oxigraph", "store_file": store_path})
    n_written = len(kg.graph)
    assert n_written > 0
    kg.close_store()

    kg2 = KnowledgeGraph(
        store="Oxigraph",
        store_file=store_path,
        identifier="http://default_graph",
    )
    n_reopen = len(kg2.graph)
    assert n_reopen == n_written, (
        f"Persistence failed: wrote {n_written} triples, reopened {n_reopen}"
    )


def test_oxigraph_purge(tmp_store_dir):
    """purge() empties the Oxigraph store."""
    store_path = os.path.join(tmp_store_dir, "kg_purge")
    kg = _build_sample_kg({"store": "Oxigraph", "store_file": store_path})
    assert len(kg.graph) > 0
    kg.purge(force=True)
    assert len(kg.graph) == 0


def test_oxigraph_multiple_structures():
    """In-memory Oxigraph KG: multiple structures can be added."""
    kg = KnowledgeGraph(store="Oxigraph")
    build.bulk("Fe", graph=kg)
    build.bulk("Cu", graph=kg)
    assert kg.n_samples == 2


def test_oxigraph_archive_roundtrip(tmp_store_dir):
    """archive() produces a valid tar.gz from an Oxigraph-backed KG."""
    kg = _build_sample_kg({"store": "Oxigraph"})
    arc_path = os.path.join(tmp_store_dir, "archive_out")
    kg.archive(arc_path)
    assert os.path.exists(arc_path + ".tar.gz")
    assert os.path.getsize(arc_path + ".tar.gz") > 0
