from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field, field_validator, PrivateAttr

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
    MDO,
    UNSAFEASMO,
)
from atomrdf.graph import KnowledgeGraph
import atomrdf.json_io as json_io
import os


class Property(DataProperty):
    graph: Optional[KnowledgeGraph] = Field(default=None)  # Private attribute for graph

    model_config = {"arbitrary_types_allowed": True}

    def _add_value(self, graph, property):
        # If the value is a list/array, serialize it to a JSON file in the
        # graph's structure_store (same approach as AtomAttribute) and add
        # identifier + path triples. Otherwise write a literal value.
        if self.value is not None:
            # handle numpy arrays and Python lists/tuples
            if isinstance(self.value, (list, tuple)) or hasattr(self.value, "tolist"):
                # ensure graph exposes a structure_store
                try:
                    store_dir = graph.structure_store
                except Exception:
                    raise RuntimeError(
                        "KnowledgeGraph has no attribute 'structure_store'"
                    )

                identifier = str(uuid.uuid4())
                datadict = {
                    identifier: {"value": self.value, "label": self.basename.lower()}
                }

                outfile = os.path.join(store_dir, str(self.id).split(":")[-1])
                json_io.write_file(outfile, datadict)
                relpath = os.path.relpath(outfile + ".json")

                graph.add(
                    (
                        property,
                        CMSO.hasIdentifier,
                        Literal(identifier, datatype=XSD.string),
                    ),
                    validate=False,
                )
                graph.add(
                    (property, CMSO.hasPath, Literal(relpath, datatype=XSD.string)),
                    validate=False,
                )
            else:
                graph.add(
                    (property, ASMO.hasValue, Literal(self.value, datatype=XSD.float)),
                    validate=False,
                )
        if self.unit is not None:
            graph.add(
                (
                    property,
                    ASMO.hasUnit,
                    URIRef(f"http://qudt.org/vocab/unit/{self.unit}"),
                ),
                validate=False,
            )

    def _create_name(self):
        name = str(uuid.uuid4())
        name = f"property:{self.basename.lower()}_{name}"
        return name

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        try:
            property = graph.create_node(
                name, getattr(ASMO, self.basename), label=self.label
            )
        except AttributeError:
            property = graph.create_node(
                name, getattr(UNSAFEASMO, self.basename), label=self.label
            )
        self._add_value(graph, property)
        return property

    @classmethod
    def from_graph(cls, graph, id):
        # get type
        typename = graph.value(id, RDF.type)
        if typename is not None:
            basename = typename.split("/")[-1]

        # get label
        label = graph.value(id, RDFS.label)
        # get value (literal) or externalised array via path+identifier
        value = graph.value(id, ASMO.hasValue)
        # get unit
        unit = graph.value(id, ASMO.hasUnit)

        # check for externalised array
        path = graph.value(id, CMSO.hasPath)
        identifier = graph.value(id, CMSO.hasIdentifier)

        cls.id = str(id)
        cls.basename = str(basename)
        cls.label = str(label) if label else None

        if value is not None:
            cls.value = float(value)
        elif path is not None and identifier is not None:
            filepath = path.toPython()
            ident = identifier.toPython()
            # open and read JSON
            with open(filepath, "r") as fin:
                data = json.load(fin)
                # defensive: ensure identifier exists
                if ident in data:
                    cls.value = data[ident]["value"]
                else:
                    raise KeyError(f"Identifier {ident} not found in {filepath}")
        else:
            cls.value = None

        cls.unit = str(unit).split("/")[-1] if unit else None
        return cls


class InputParameter(Property):
    basename: Optional[str] = Field(
        default="InputParameter", description="Basename of the property"
    )


class OutputParameter(Property):
    basename: Optional[str] = Field(
        default="OutputParameter", description="Basename of the property"
    )
    associate_to_sample: Optional[List[str]] = Field(
        default=None, description="List of sample IDs to associate the property with"
    )


class CalculatedProperty(Property):
    basename: Optional[str] = Field(
        default="CalculatedProperty", description="Basename of the property"
    )
    associate_to_sample: Optional[List[str]] = Field(
        default=None, description="List of sample IDs to associate the property with"
    )
