import pytest
from atomrdf import KnowledgeGraph, System

def test_structuregraph():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	assert(sys.sample != None)

	sys = System.create.lattice.bcc(element="Fe", graph=s)
	assert(sys.sample != None)

	sys = System.read.file("tests/al_data/Al.poscar", format="poscar", graph=s)
	assert(sys.sample != None)

	sys = System.create.defect.grain_boundary(axis=[0,0,1], 
                        sigma=5, 
                        gb_plane=[3, -1, 0],
                        element='Fe',
                        graph=s)

	assert(sys.sample != None)




