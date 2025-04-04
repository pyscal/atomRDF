import pytest
import os
from atomrdf import KnowledgeGraph, System
from atomrdf.namespace import CMSO, PLDO
import shutil
import atomrdf.build as build

def test_sqlalchemy():
	s = KnowledgeGraph(store='SQLAlchemy', store_file='aa')
	sys = build.bulk("Fe", graph=s)