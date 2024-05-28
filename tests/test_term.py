import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil

def test_lt():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = term < 2
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue<"2"^^xsd:int)'

def test_eq():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = (term == 2)
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue="2"^^xsd:int)'	

def test_gt():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = term > 2
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue>"2"^^xsd:int)'

def test_lte():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = term <= 2
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue<="2"^^xsd:int)'

def test_gte():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = term >= 2
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue>="2"^^xsd:int)'

def test_ne():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = kg.ontology.terms.cmso.hasNumberOfAtoms
	term = term != 2
	assert term._condition == '(?cmso:hasNumberOfAtomsvalue!="2"^^xsd:int)'

def test_and():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = (kg.ontology.terms.cmso.hasVolume > 2) & (kg.ontology.terms.cmso.hasVolume < 4)
	assert term._condition == '((?cmso:hasVolumevalue>"2"^^xsd:float)&&(?cmso:hasVolumevalue<"4"^^xsd:float))'

def test_or():
	kg = KnowledgeGraph()
	struct_Fe = System.create.element.Fe(graph=kg)

	term = (kg.ontology.terms.cmso.hasVolume > 2) | (kg.ontology.terms.cmso.hasVolume < 4)
	assert term._condition == '((?cmso:hasVolumevalue>"2"^^xsd:float)||(?cmso:hasVolumevalue<"4"^^xsd:float))'