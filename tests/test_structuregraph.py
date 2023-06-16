import pytest
from pyscal_rdf import StructureGraph

def test_structuregraph():
	s = StructureGraph()
	sys = s.create_element("Fe")
	assert(s.sample != None)

	sys = s.create_structure("bcc", element="Fe")
	assert(s.sample != None)

	sys = s.read_structure("tests/al_data/Al.poscar", format="poscar")
	assert(s.sample != None)

	sys = s.create_grain_boundary(axis=[0,0,1], 
                        sigma=5, 
                        gb_plane=[3, -1, 0],
                        element='Fe')
	assert(s.sample != None)




