import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil


def test_network_draw():
	s = KnowledgeGraph()
	sys = System.create.element.Fe(graph=s)
	assert s.ontology.draw() != None