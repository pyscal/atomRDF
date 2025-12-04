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


class Method(BaseModel, TemplateMixin):
    basename: str
    pid: Optional[str] = Field(default=None, description="PID of the method")


class MolecularStatics(Method):
    basename: str = "MolecularStatics"
    pid: str = ASMO.MolecularStatics.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.MolecularStatics)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class MolecularDynamics(Method):
    basename: str = "MolecularDynamics"
    pid: str = ASMO.MolecularDynamics.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.MolecularDynamics)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class DensityFunctionalTheory(Method):
    basename: str = "DensityFunctionalTheory"
    pid: str = ASMO.DensityFunctionalTheory.uri

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_method"
        method = graph.create_node(main_id, ASMO.DensityFunctionalTheory)
        return method

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


method_map = {
    "MolecularDynamics": MolecularDynamics,
    "MolecularStatics": MolecularStatics,
    "DensityFunctionalTheory": DensityFunctionalTheory,
}
