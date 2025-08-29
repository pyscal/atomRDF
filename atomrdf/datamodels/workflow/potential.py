from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field, field_validator
from atomrdf.datamodels.basemodels import (
    TemplateMixin,
    DataProperty,
    RDFMixin,
    BaseModel,
    Activity,
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


class InteratomicPotential(BaseModel, TemplateMixin):
    pid: str = ASMO.InteratomicPotential.uri
    uri: Optional[str] = Field(default=None, description="URI of the potential")
    potential_type: str = Field(
        default="InteratomicPotential", description="Type of the potential"
    )

    def _add_potential(self, potential, graph):
        if self.uri:
            graph.add(
                (potential, CMSO.hasReference, Literal(self.uri, datatype=XSD.string))
            )
        if self.label:
            graph.add((potential, RDFS.label, Literal(self.label, datatype=XSD.string)))

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_potential"
        potential = graph.create_node(main_id, ASMO.InteratomicPotential)
        self._add_potential(potential)
        return potential


class ModifiedEmbeddedAtomModel(InteratomicPotential):
    pid: str = ASMO.ModifiedEmbeddedAtomModel.uri
    potential_type: str = Field(
        default="ModifiedEmbeddedAtomModel", description="Type of the potential"
    )

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_potential"
        potential = graph.create_node(main_id, ASMO.ModifiedEmbeddedAtomModel)
        self._add_potential(potential)
        return potential


class EmbeddedAtomModel(InteratomicPotential):
    pid: str = str(ASMO.EmbeddedAtomModel.uri)
    potential_type: str = Field(
        default="EmbeddedAtomModel", description="Type of the potential"
    )

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_potential"
        potential = graph.create_node(main_id, ASMO.EmbeddedAtomModel)
        self._add_potential(potential)
        return potential


class LennardJonesPotential(InteratomicPotential):
    pid: str = ASMO.LennardJonesPotential.uri
    potential_type: str = Field(
        default="LennardJonesPotential", description="Type of the potential"
    )

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_potential"
        potential = graph.create_node(main_id, ASMO.LennardJonesPotential)
        self._add_potential(potential)
        return potential


class MachineLearningPotential(InteratomicPotential):
    pid: str = ASMO.MachineLearningPotential.uri
    potential_type: str = Field(
        default="MachineLearningPotential", description="Type of the potential"
    )

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_potential"
        potential = graph.create_node(main_id, ASMO.MachineLearningPotential)
        self._add_potential(potential)
        return potential
