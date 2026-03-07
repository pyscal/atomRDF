from typing import ClassVar, List, Optional, Union
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
    pid: str = str(ASMO.MolecularStatics)
    base_uriref: ClassVar[URIRef] = URIRef("method_MolecularStatics")

    def to_graph(self, graph):
        if graph.persistent_members["MolecularStatics"] is not None:
            return graph.persistent_members["MolecularStatics"]
        else:
            graph.persistent_members["MolecularStatics"] = self.base_uriref
        return self.base_uriref

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class MolecularDynamics(Method):
    basename: str = "MolecularDynamics"
    pid: str = str(ASMO.MolecularDynamics)
    base_uriref: ClassVar[URIRef] = URIRef("method_MolecularDynamics")

    def to_graph(self, graph):
        if graph.persistent_members["MolecularDynamics"] is not None:
            return graph.persistent_members["MolecularDynamics"]
        else:
            graph.persistent_members["MolecularDynamics"] = self.base_uriref
        return self.base_uriref

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


class DensityFunctionalTheory(Method):
    basename: str = "DensityFunctionalTheory"
    pid: str = str(ASMO.DensityFunctionalTheory)
    base_uriref: ClassVar[URIRef] = URIRef("method_DensityFunctionalTheory")

    def to_graph(self, graph, main_id):
        if graph.persistent_members["DensityFunctionalTheory"] is not None:
            return graph.persistent_members["DensityFunctionalTheory"]
        else:
            graph.persistent_members["DensityFunctionalTheory"] = self.base_uriref
        return self.base_uriref

    @classmethod
    def from_graph(cls, graph, id):
        label = graph.get_label(id)
        return cls(label=label)


method_map = {
    "MolecularDynamics": MolecularDynamics,
    "MolecularStatics": MolecularStatics,
    "DensityFunctionalTheory": DensityFunctionalTheory,
}
