import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO

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

def test_visualise():
	s = KnowledgeGraph()
	sys = System.create.element.Cr(graph=s)

	styledict = {"edgecolor": "#D9D9D9",
	"BNode": {"color": "#263238"}}
	vis = s.visualise(styledict=styledict)
	assert(vis != None)

def test_logger():
	if os.path.exists('tests/atomrdf.log'):
		os.remove('tests/atomrdf.log')
	s = KnowledgeGraph(enable_log=True)
	s.log('testing-logger')
	assert str(type(s.log).__name__) == "method"
 
def test_add_structure():
	s = KnowledgeGraph()
	sys = System.create.element.Fe()
	s.add_structure(sys)
	assert sys.sample in s.samples

def test_add_cross_triple():
	
	s = KnowledgeGraph(enable_log=True)
	sys = System.create.element.Fe(graph=s)
	status, _ = s._check_domain_if_uriref((sys.material, CMSO.hasDefect, PLDO.AntiphaseBoundary))
	assert status == True
	

def test_add_quantity():
	s = KnowledgeGraph(enable_log=True)
	sys = System.create.element.Fe(graph=s)
	s.add_calculated_quantity(sys.sample,
		'Energy',
		str(23),
		unit='eV')
	cp = s.value(sys.sample, CMSO.hasCalculatedProperty)
	val = s.value(cp, CMSO.hasValue)
	assert val.toPython() == '23'

#def test_add_domain_ontoterm():
#	from atomrdf.namespace import CMSO, PLDO
#	s = KnowledgeGraph()
#	sys = System.create.element.Fe(graph=s)
#	status, _ = s._check_domain_if_ontoterm((CMSO.Material, CMSO.hasDefect, PLDO.AntiphaseBoundary))
#	assert status == True



