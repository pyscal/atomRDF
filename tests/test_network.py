import pytest
import os
from atomrdf import KnowledgeGraph

import shutil
import atomrdf.build as build


def test_network_draw():
    s = KnowledgeGraph()
    sys = build.bulk("Fe", graph=s)
    if s.ontology is None:
        pytest.skip("Ontology not available (network/purls unavailable)")
    assert s.ontology.draw() != None
