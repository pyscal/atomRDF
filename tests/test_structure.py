import pytest
import os
import numpy as np
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil

def test_custom():
	kg = KnowledgeGraph()
	struct = System.create.lattice.custom([[0.   , 0.   , 0.   ],
	       [1.435, 1.435, 1.435]],
	    [1, 1], 
	    [[2.87, 0.  , 0.  ],
	       [0.  , 2.87, 0.  ],
	       [0.  , 0.  , 2.87]],
	    lattice_constant = 2.87, 
	    element='Fe',
	    graph=kg)
	assert kg.value(struct.sample, CMSO.hasNumberOfAtoms).toPython() == 2

def test_dislocation():
	slip_direction = np.array([1, 0, -1])
	slip_plane = np.array([1, 1, 1])
	slip_system = [slip_direction, slip_plane]
	burgers_vector = 0.5
	dislocation_line = np.array([1, 0, -1])
	elastic_constant_dict = {'C11': 169, 'C12': 122, 'C44': 75.4}
	kg = KnowledgeGraph()
	sys = System.create.defect.dislocation(slip_system,
                                        dislocation_line,
                                        elastic_constant_dict,
                                        burgers_vector=burgers_vector,
                                        element='Cu',
                                        dislocation_type='monopole',
                                        graph=kg,
                                        ) 
	assert sys.natoms == 96

	sys = System.create.defect.dislocation(slip_system,
                                        dislocation_line,
                                        elastic_constant_dict,
                                        burgers_vector=burgers_vector,
                                        element='Cu',
                                        dislocation_type='periodicarray',
                                        graph=kg,
                                        ) 
	assert sys.natoms == 96

	res = kg.query_sample(kg.ontology.terms.ldo.ScrewDislocation)
	assert res is not None


def test_read_in():
	kg = KnowledgeGraph()
	struct = System.read.file('tests/conf.dump', graph=kg, lattice='bcc', lattice_constant=2.861)
	assert kg.n_samples == 1

def test_delete():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	sys.delete(indices=[0])
	assert sys.natoms == 1
	ss, n= s.get_sample(sys.sample, no_atoms=True)
	assert n==1

	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	del sys[0]
	assert sys.natoms == 1
	ss, n= s.get_sample(sys.sample, no_atoms=True)
	assert n==1

def test_substitute():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	sys.substitute_atoms('Li', indices=[0])
	species = s.value(sys.sample, CMSO.hasSpecies)
	elements = [k[2] for k in s.triples((species, CMSO.hasElement, None))]
	assert len(elements) == 2

def test_interstitials():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	sys = sys.add_interstitial_impurities(['Li', 'Au'], void_type='tetrahedral')
	species = s.value(sys.sample, CMSO.hasSpecies)
	elements = [k[2] for k in s.triples((species, CMSO.hasElement, None))]
	assert len(elements) == 3

	sys = System.create.element.Fe(graph=s)
	sys = sys.add_interstitial_impurities(['Li', 'Au'], void_type='octahedral')
	species = s.value(sys.sample, CMSO.hasSpecies)
	elements = [k[2] for k in s.triples((species, CMSO.hasElement, None))]
	assert len(elements) == 3

def test_gb():
	kg = KnowledgeGraph()
	struct_gb_1 = System.create.defect.grain_boundary(axis=[0,0,1], 
                        sigma=5, 
                        gb_plane=[3, -1, 0],
                        element='Fe',
                        graph=kg)
	res = kg.query_sample(kg.ontology.terms.pldo.GrainBoundary)
	assert len(res.AtomicScaleSample.values) == 1
	
	new = struct_gb_1.repeat((2,2,2))
	res = kg.query_sample(kg.ontology.terms.pldo.GrainBoundary)
	assert len(res.AtomicScaleSample.values) == 2

	ss = kg.get_sample(new.sample)
	res = ss.query_sample(ss.ontology.terms.pldo.GrainBoundary)
	assert len(res.AtomicScaleSample.values) == 1

def test_sf():
	kg = KnowledgeGraph()
	slip_plane = np.array([0,0,0,1])
	slip_direction = np.array([1,0,-1,0])
	sf = System.create.defect.stacking_fault(slip_plane, 
    1.0,
    slip_direction_a=slip_direction,
    element="Mg", 
    repetitions=(2,2,2),
    vacuum=10.0,
    graph=kg)
	assert sf.natoms == 96

	slip_plane = np.array([1,1,1])
	sf, sfa, ssa, fsa = System.create.defect.stacking_fault(slip_plane, 
		0.5,
		element="Cu", 
		repetitions=(1,1,1),
		vacuum=0.0,
		graph=kg,
		return_atomman_dislocation=True)
	assert sfa.system.natoms == 48
	assert sfa.abovefault.sum() == 24
	assert np.abs(sfa.faultpos_rel - 0.5) < 1e-6
	assert np.abs(sfa.faultpos_cart - 12.505406830647294) < 1e-6