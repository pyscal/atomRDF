import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil


def test_sqlalchemy():
	s = KnowledgeGraph(store='SQLAlchemy', store_file='aa')
	sys = System.create.element.Fe(graph=s)