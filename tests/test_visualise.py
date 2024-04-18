import pytest
import os
import numpy as np
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil
import atomrdf.visualize as viz


def test_switch_box():
	assert viz._switch_box("box") == "rectangle"

def test_fix_id():
	assert viz._fix_id('hello', 'random') == 'hello'
	assert viz._fix_id('hello', 'Literal') != 'hello'

def test_visualise():
	s = KnowledgeGraph()
	sys = System.create.element.Cr(graph=s)

	styledict = {"edgecolor": "#D9D9D9",
	"BNode": {"color": "#263238"}}
	vis = s.visualise(styledict=styledict)
	assert(vis != None)

	s = KnowledgeGraph()
	sys = System.create.element.Cr(graph=s)
	vis = s.visualise(workflow_view=True, hide_types=True)
	assert(vis != None)
