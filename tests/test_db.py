from pyscal_rdf import StructureGraph
import pytest


def test_db():
	g = StructureGraph(store="SQLAlchemy", store_file="testfile.db")
	struct_Fe = g.create_element("Fe")
	assert g.n_samples == 1

	k = StructureGraph(store="SQLAlchemy", store_file="testfile.db")
	assert k.n_samples == 1
