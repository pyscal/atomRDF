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
    MDO,
)


class Property(DataProperty):
    def _add_value(self, graph, property):
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


    def _create_name(self):
        name = str(uuid.uuid4())
        name = f"property:{self.basename.lower()}_{name}"
        return name

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, getattr(ASMO, self.basename), label=self.label)
        self._add_value(graph, property)
        return property
        
    @classmethod
    def from_graph(cls, graph, id):
        #get type
        typename = graph.value(id, RDF.type)
        if typename is not None:
            basename = typename.split("/")[-1]
        
        # get label
        label = graph.value(id, RDFS.label)
        # get value
        value = graph.value(id, ASMO.hasValue)
        # get unit
        unit = graph.value(id, ASMO.hasUnit)

        cls.id = str(id)
        cls.basename = str(basename)
        cls.label = str(label) if label else None
        cls.value = float(value) if value else None
        cls.unit = str(unit).split("/")[-1] if unit else None
        return cls

class InputParameter(Property):
    basename: Optional[str] = Field(
        default='InputParameter', description="Basename of the property"
    )    

class OutputParameter(Property):
    basename: Optional[str] = Field(
        default='OutputParameter', description="Basename of the property"
    )    
    associate_to_sample: Optional[bool] = Field(
        default=True, description="Whether to associate the property to the sample"
    )
    
class CalculatedProperty(Property):
    basename: Optional[str] = Field(
        default='CalculatedProperty', description="Basename of the property"
    )
    associate_to_sample: Optional[bool] = Field(
        default=True, description="Whether to associate the property to the sample"
    )    
    