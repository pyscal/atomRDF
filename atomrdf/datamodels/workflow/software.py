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


class SoftwareAgent(BaseModel, TemplateMixin):
    uri: Optional[str] = Field(default=None, description="URI of the software agent")
    version: Optional[str] = Field(
        default=None, description="Version of the software agent"
    )

    def to_graph(
        self,
        graph,
    ):
        agent = graph.create_node(self.uri, PROV.SoftwareAgent, label=self.label)
        return agent
        # if method:
        #    graph.add((method, PROV.wasAssociatedWith, agent))
        # if workflow_agent:
        #    graph.add((workflow_agent, PROV.actedOnBehalfOf, agent))

    @classmethod
    def from_graph(cls, graph, id):
        uri = graph.value(id, RDFS.type)
        label = graph.value(id, RDFS.label)
        cls.uri = str(uri) if uri else None
        cls.label = str(label) if label else None
        return cls
