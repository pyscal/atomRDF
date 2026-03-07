"""
Tests for Dataset, Publication, and Creator datamodels.
"""

import pytest
from rdflib import URIRef, Literal
from atomrdf import KnowledgeGraph
from atomrdf.datamodels.dataset import Creator, Dataset, Publication
from atomrdf.namespace import DCAT, DCTERMS, FOAF


@pytest.fixture
def kg():
    return KnowledgeGraph()


# ---------------------------------------------------------------------------
# Creator
# ---------------------------------------------------------------------------


def test_creator_to_graph_with_id(kg):
    creator = Creator(id="https://orcid.org/0000-0001-2345-6789", name="Alice Smith")
    uri = creator.to_graph(kg)
    assert uri == URIRef("https://orcid.org/0000-0001-2345-6789")
    assert (uri, FOAF.name, Literal("Alice Smith")) in kg.graph


def test_creator_to_graph_without_id(kg):
    creator = Creator(name="Bob Jones")
    uri = creator.to_graph(kg)
    assert str(uri).startswith("person:")
    assert (uri, FOAF.name, Literal("Bob Jones")) in kg.graph


def test_creator_from_graph(kg):
    creator = Creator(id="https://orcid.org/0000-0000-0000-0001", name="Carol White")
    uri = creator.to_graph(kg)
    recovered = Creator.from_graph(kg, uri)
    assert recovered.name == "Carol White"
    assert recovered.id == "https://orcid.org/0000-0000-0000-0001"


# ---------------------------------------------------------------------------
# Publication
# ---------------------------------------------------------------------------


def test_publication_to_graph_with_id(kg):
    pub = Publication(
        id="https://doi.org/10.1016/j.actamat.2024.00001",
        identifier="10.1016/j.actamat.2024.00001",
        title="Grain boundary energies",
    )
    uri = pub.to_graph(kg)
    assert uri == URIRef("https://doi.org/10.1016/j.actamat.2024.00001")
    assert (
        uri,
        DCTERMS.identifier,
        Literal("10.1016/j.actamat.2024.00001"),
    ) in kg.graph
    assert (uri, DCTERMS.title, Literal("Grain boundary energies")) in kg.graph


def test_publication_to_graph_without_id(kg):
    pub = Publication(identifier="10.1016/j.actamat.2024.99999")
    uri = pub.to_graph(kg)
    assert str(uri).startswith("publication:")
    assert (
        uri,
        DCTERMS.identifier,
        Literal("10.1016/j.actamat.2024.99999"),
    ) in kg.graph


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


def test_dataset_to_graph(kg):
    ds = Dataset(
        identifier="https://doi.org/10.5281/zenodo.1234567",
        title="Test dataset",
        creators=[Creator(id="https://orcid.org/0000-0001-0001-0001", name="Dave K")],
        publication=Publication(
            id="https://doi.org/10.1234/paper.001",
            identifier="10.1234/paper.001",
            title="A paper",
        ),
    )
    uri = ds.to_graph(kg)
    assert uri == URIRef("https://doi.org/10.5281/zenodo.1234567")
    assert (uri, DCAT.Dataset, None) not in kg.graph  # type triple uses RDF.type
    assert (uri, DCTERMS.title, Literal("Test dataset")) in kg.graph
    assert (
        uri,
        DCTERMS.identifier,
        Literal("https://doi.org/10.5281/zenodo.1234567"),
    ) in kg.graph
    # Creator link
    person_uri = URIRef("https://orcid.org/0000-0001-0001-0001")
    assert (uri, DCTERMS.creator, person_uri) in kg.graph
    # Publication link
    pub_uri = URIRef("https://doi.org/10.1234/paper.001")
    assert (uri, DCTERMS.isReferencedBy, pub_uri) in kg.graph


def test_dataset_samples_isPartOf(kg):
    sample_uri = "https://example.org/sample/001"
    ds = Dataset(
        identifier="https://doi.org/10.5281/zenodo.9999999",
        samples=[sample_uri],
    )
    uri = ds.to_graph(kg)
    assert (URIRef(sample_uri), DCTERMS.isPartOf, uri) in kg.graph


def test_dataset_from_dict(kg):
    data = {
        "identifier": "https://doi.org/10.5281/zenodo.7654321",
        "title": "From dict dataset",
        "creators": [{"id": "https://orcid.org/0000-0002-0002-0002", "name": "Eve A"}],
        "publication": {
            "id": "https://doi.org/10.9999/paper.002",
            "identifier": "10.9999/paper.002",
            "title": "Another paper",
        },
        "samples": [],
    }
    ds = Dataset.from_dict(data)
    assert ds.identifier == "https://doi.org/10.5281/zenodo.7654321"
    assert ds.title == "From dict dataset"
    assert len(ds.creators) == 1
    assert ds.creators[0].name == "Eve A"
    assert ds.publication is not None
    assert ds.publication.identifier == "10.9999/paper.002"
    uri = ds.to_graph(kg)
    assert (uri, DCTERMS.title, Literal("From dict dataset")) in kg.graph


# ---------------------------------------------------------------------------
# WorkflowParser integration
# ---------------------------------------------------------------------------


def test_workflow_parser_dataset_key():
    """WorkflowParser should handle a top-level 'dataset' key."""
    from atomrdf.io.workflow_parser import WorkflowParser

    data = {
        "dataset": {
            "identifier": "https://doi.org/10.5281/zenodo.1111111",
            "title": "Parser integration test",
            "creators": [{"name": "Frank B"}],
        }
    }
    parser = WorkflowParser()
    result = parser.parse(data)
    assert result["dataset_uri"] == "https://doi.org/10.5281/zenodo.1111111"
    ds_uri = URIRef("https://doi.org/10.5281/zenodo.1111111")
    assert (
        ds_uri,
        DCTERMS.title,
        Literal("Parser integration test"),
    ) in parser.kg.graph


def test_workflow_parser_dataset_sample_resolution():
    """Sample IDs in dataset.samples should be resolved via the sample_map."""
    from atomrdf.io.workflow_parser import WorkflowParser

    # We test the resolution logic without adding a real structure — just inject
    # a pre-populated sample_map via the parser.
    parser = WorkflowParser()
    parser.sample_map = {"sample_A": "https://example.org/sample/real_uri"}

    data = {
        "dataset": {
            "identifier": "https://doi.org/10.5281/zenodo.2222222",
            "samples": ["sample_A"],
        }
    }
    result = parser.parse(data)
    ds_uri = URIRef("https://doi.org/10.5281/zenodo.2222222")
    assert (
        URIRef("https://example.org/sample/real_uri"),
        DCTERMS.isPartOf,
        ds_uri,
    ) in parser.kg.graph
