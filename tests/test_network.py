import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil
import atomrdf.build as build

def test_network_draw():
	s = KnowledgeGraph()
	sys = build.bulk("Fe", graph=s)
	assert s.ontology.draw() != None