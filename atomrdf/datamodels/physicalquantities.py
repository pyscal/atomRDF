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


class PhysicalQuantity(DataProperty, TemplateMixin):
    """
    A class to represent a physical quantity with its value, unit, and associated metadata.
    """

    def to_graph(
        self,
        graph,
    ):
        # this knows just the serialisation of itself.
        name = str(uuid.uuid4())
        if self.basename:
            name = f"{self.basename.lower()}:{name}"
        else:
            name = f"property:{name}"
        self.id = name
        property = graph.create_node(
            name, getattr(ASMO, self.basename), label=self.label
        )
        if self.value is not None:
            graph.add(
                (property, ASMO.hasValue, Literal(self.value, datatype=XSD.float))
            )
        if self.unit is not None:
            graph.add(
                (
                    property,
                    ASMO.hasUnit,
                    URIRef(f"http://qudt.org/vocab/unit/{self.unit}"),
                )
            )

    @classmethod
    def from_graph(cls, graph, id):
        # recover PID
        pid = graph.value(id, RDFS.type)
        pid = str(pid) if pid else None
        basename = str(pid).split("/")[-1] if pid else None

        # get label
        label = graph.value(id, RDFS.label)
        # get value
        value = graph.value(id, ASMO.hasValue)
        # get unit
        unit = graph.value(id, ASMO.hasUnit)
        cls.id = str(id)
        cls.pid = str(pid) if pid else None
        cls.basename = str(basename) if basename else None
        cls.label = str(label) if label else None
        cls.value = float(value) if value else None
        cls.unit = str(unit).split("/")[-1] if unit else None
        return cls
