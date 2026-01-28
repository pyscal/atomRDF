from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field
from atomrdf.datamodels.basemodels import (
    TemplateMixin,
    DataProperty,
    RDFMixin,
    BaseModel,
)
from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.namespace import (
    CMSO,
    LDO,
    PLDO,
    PODO,
    CDCO,
    PROV,
    Literal,
    ASMO,
)


class DegreeOfFreedom(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the degree of freedom")


class AtomicPositionRelaxation(DegreeOfFreedom):
    pid: str = ASMO.AtomicPositionRelaxation.uri
    basename: str = "AtomicPositionRelaxation"

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return ASMO.AtomicPositionRelaxation


class CellVolumeRelaxation(DegreeOfFreedom):
    pid: str = ASMO.CellVolumeRelaxation.uri
    basename: str = "CellVolumeRelaxation"

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return ASMO.CellVolumeRelaxation


class CellShapeRelaxation(DegreeOfFreedom):
    pid: str = ASMO.CellShapeRelaxation.uri
    basename: str = "CellShapeRelaxation"

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return ASMO.CellShapeRelaxation


dof_map = {
    "AtomicPositionRelaxation": AtomicPositionRelaxation,
    "CellVolumeRelaxation": CellVolumeRelaxation,
    "CellShapeRelaxation": CellShapeRelaxation,
}
