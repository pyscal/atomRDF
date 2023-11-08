import pytest
from pyscal_rdf import StructureGraph

def test_structuregraph():
	s = StructureGraph()
	sys = s.create.element.Fe()
	assert(s.sample != None)

	sys = s.create.lattice.bcc(element="Fe")
	assert(s.sample != None)

	sys = s.read_structure("tests/al_data/Al.poscar", format="poscar")
	assert(s.sample != None)

	sys = s.create.defect.grain_boundary(axis=[0,0,1], 
                        sigma=5, 
                        gb_plane=[3, -1, 0],
                        element='Fe')
	assert(s.sample != None)




