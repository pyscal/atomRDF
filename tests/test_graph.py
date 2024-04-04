import pytest
from pyscal_rdf import KnowledgeGraph, System

def test_structuregraph():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)

	vis = s.visualise()
	assert(vis != None)

	vis = s.visualise()
	assert(vis != None)

	s.write("temp.ttl", format="turtle")
	s = KnowledgeGraph(graph_file="temp.ttl")

	sys = System.create.element.Fe(graph=s)
	sys.add_vacancy(0.5, number=1)

	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	assert s.n_samples == 1
	#res = s.query_sample("NumberOfAtoms", 2)
	#assert(len(res) == 1)