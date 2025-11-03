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


class Algorithm(BaseModel, TemplateMixin):
    basename: str
    pid: Optional[str] = Field(default=None, description="PID of the algorithm")
    label: Optional[str] = Field(default=None, description="Label of the algorithm")


class EquationOfStateFit(Algorithm):
    basename: str = "EquationOfStateFit"
    pid: str = ASMO.EquationOfStateFit.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.EquationOfStateFit)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class QuasiHarmonicApproximation(Algorithm):
    basename: str = "QuasiHarmonicApproximation"
    pid: str = ASMO.QuasiHarmonicApproximation.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.QuasiHarmonicApproximation)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class ThermodynamicIntegration(Algorithm):
    basename: str = "ThermodynamicIntegration"
    pid: str = ASMO.ThermodynamicIntegration.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.ThermodynamicIntegration)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class ANNNIModel(Algorithm):
    basename: str = "ANNNIModel"
    pid: str = ASMO.ANNImodel.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.ANNImodel)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


algorithm_map = {
    "EquationOfStateFit": EquationOfStateFit,
    "QuasiHarmonicApproximation": QuasiHarmonicApproximation,
    "ThermodynamicIntegration": ThermodynamicIntegration,
    "ANNIModel": ANNNIModel,
}
