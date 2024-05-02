import pytest
import os
import numpy as np
from atomrdf import KnowledgeGraph, System, Workflow
from atomrdf.namespace import CMSO, PLDO
import shutil

def test_wf_creation():
	kg = KnowledgeGraph()
	#kg.enable_workflow(pr, workflow_environment='pyiron')

