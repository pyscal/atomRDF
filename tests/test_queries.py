import pytest
from pyscal_rdf import StructureGraph
from pyscal_rdf.queries import Query

def test_query():
	s = StructureGraph()
	sys = s.create_element("Fe")
	q = Query()

	assert(len(q.sparql.sample_by_latticesystem(s, "Q851536")) == 1)

	struct_gb_1 = s.create_grain_boundary(axis=[0,0,1], 
                        sigma=5, 
                        gb_plane=[3, -1, 0],
                        element='Fe')
	assert(len(q.sparql.sample_by_defect(s, "symmetric tilt")) == 1)
	assert(len(q.sparql.sample_by_sigma(s, 5)) == 1)
	assert(len(q.python.sample_by_altname(s, "bcc")) == 2)




