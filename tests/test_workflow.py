import pytest
import os
import numpy as np
from atomrdf import KnowledgeGraph, System, Workflow
from atomrdf.namespace import CMSO, PLDO
import shutil

def test_wf_creation():
	s = KnowledgeGraph()
	wf = Workflow(s)

def test_lattice_props():
	s = KnowledgeGraph()
	wf = Workflow(s)
	parent_sys = System.create.element.Cr(graph=s)
	#add some defects
	parent_sys.delete(indices=[0])
	child_sys = System.create.element.Cr(graph=s)
	#trick workflow and add info
	wf.parent_sample = parent_sys.sample
	wf.structure = child_sys
	wf.sample = child_sys.sample
	wf.add_structural_relation()


def test_add_method():
	for pot in ['meam', 'eam', 'lj', 'ace', 'nequip']:
		s = KnowledgeGraph()
		wf = Workflow(s)
		parent_sys = System.create.element.Cr(graph=s)
		#add some defects
		parent_sys.delete(indices=[0])
		child_sys = System.create.element.Cr(graph=s)
		#trick workflow and add info
		wf.parent_sample = parent_sys.sample
		wf.structure = child_sys
		wf.sample = child_sys.sample
		wf.add_structural_relation()
		wf.main_id = 2314
		wf.mdict = {"method": "MolecularDynamics",
		"temperature": 100,
		"pressure": 0,
		"dof": ["AtomicPosition", "CellVolume"],
		"ensemble": "IsothermalisobaricEnsemble",
		"id": 2314,
		"potential": {"uri": "https://doi.org/xxx",
		  "type": pot,
		  "label": "string" },
		"workflow_manager": {"uri": "xxxx",
		  "label": "pyiron"},
		"software": [ {"uri": "xxxx", "label": "lammps"},],
		"outputs": [{"label": "TotalEnergy", "value": 2.301, "unit": "EV", 
		"associate_to_sample": True}],
		"inputs": [ {"label": "AnotherInput", "value": 0.1, "unit": None },]
		}
		wf.add_method()
