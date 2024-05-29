import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO, ASMO
import shutil

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
	assert sys.sample in s.sample_ids

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
	val = s.value(cp, ASMO.hasValue)
	assert val.toPython() == '23'

	insp = s.inspect_sample(sys.sample)
	assert 'Im-3m' in insp
	assert '23' in insp

def test_archive():
	s = KnowledgeGraph(enable_log=True)
	sys = System.create.element.Fe(graph=s)
	sys = System.create.element.Cu(graph=s)
	if os.path.exists('test_archive.tar.gz'):
		os.remove('test_archive.tar.gz')
	if os.path.exists('test_archive'):
		shutil.rmtree('test_archive')
	s.archive('test_archive')
	assert os.path.exists('test_archive.tar.gz')

	s = KnowledgeGraph.unarchive('test_archive.tar.gz')
	assert s.n_samples == 2
	os.remove('test_archive.tar.gz')
	shutil.rmtree('test_archive')

def test_sparql_query():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)
	struct_Si = System.create.element.Si(graph=kg)
	struct_l12 = System.create.lattice.l12(element=['Al', 'Ni'], 
                         lattice_constant=3.57, graph=kg)
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
	res = kg.query(query)
	assert res.symbol.values[0].toPython() == 'Pm-3m'

	res = kg.query_sample(kg.ontology.terms.cmso.hasAltName=='bcc')
	assert res.hasAltNamevalue.values[0].toPython() == 'bcc'

	res = kg.query_sample(kg.ontology.terms.cmso.hasAltName=='bcc',
	             enforce_types=True)
	assert res.hasAltNamevalue.values[0].toPython() == 'bcc'

def test_extract_sample():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)
	sample_graph, no_atoms = kg.get_sample(struct_Fe.sample, no_atoms=True)
	assert no_atoms == 2
	assert sample_graph.sample_ids[0] == struct_Fe.sample

	struct = kg.get_system_from_sample(struct_Fe.sample)
	assert len(struct.atoms.positions) == 2
	assert struct.graph is not None

	kg.to_file(struct_Fe.sample, filename='POSCAR')
	assert os.path.exists('POSCAR')
	os.remove('POSCAR')

	kg.to_file(struct_Fe.sample, filename='POSCAR', format='cif')
	assert os.path.exists('POSCAR')
	os.remove('POSCAR')

#def test_add_domain_ontoterm():
#	from atomrdf.namespace import CMSO, PLDO
#	s = KnowledgeGraph()
#	sys = System.create.element.Fe(graph=s)
#	status, _ = s._check_domain_if_ontoterm((CMSO.Material, CMSO.hasDefect, PLDO.AntiphaseBoundary))
#	assert status == True



