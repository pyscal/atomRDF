import pytest
from atomrdf import KnowledgeGraph, System
import atomrdf.build as build
import atomrdf.io as aio


def test_structuregraph():
    s = KnowledgeGraph()
    sys = build.bulk("Fe", graph=s)
    assert sys.sample != None

    sys = build.bulk("Fe", structure="bcc", graph=s)
    assert sys.sample != None

    sys = aio.read("tests/al_data/Al.poscar", format="poscar", graph=s)
    assert sys.sample != None

    sys = build.defect.grain_boundary(
        axis=[0, 0, 1],
        sigma=5,
        gb_plane=[3, -1, 0],
        element="Fe",
        graph=s,
        backend="inbuilt",
    )

    assert sys.sample != None

    sys = build.defect.grain_boundary(
        axis=[0, 0, 1],
        sigma=5,
        gb_plane=[3, -1, 0],
        element="Fe",
        graph=s,
        backend="aimsgb",
    )

    assert sys.sample != None
