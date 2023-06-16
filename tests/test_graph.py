import pytest
from pyscal_rdf import StructureGraph

def test_structuregraph():
	s = StructureGraph()
	sys = s.create_element("Fe")

	vis = s.visualise()
	assert(vis != None)

	vis = s.visualise(styledict={"BNode":{"fontsize":10}})
	assert(vis != None)

	s.write("temp.ttl", format="turtle")
	s = StructureGraph(graph_file="temp.ttl")

	sys = s.create_element("Fe")
	s.add_vacancy(0.5, number=1)

	s = StructureGraph()
	sys = s.create_element("Fe")
	res = s.query_sample("NumberOfAtoms", 2)
	assert(len(res) == 1)