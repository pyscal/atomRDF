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
    def _create_name(self):
        name = str(uuid.uuid4())
        basename = self.__class__.__name__
        name = f"{basename.lower()}:{name}"
        return name

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

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.InputParameter, label=self.label)
        self._add_value(graph, property)

    @classmethod
    def _get_value(cls, graph, id):
        # get label
        label = graph.value(id, RDFS.label)
        # get value
        value = graph.value(id, ASMO.hasValue)
        # get unit
        unit = graph.value(id, ASMO.hasUnit)
        cls.id = str(id)
        cls.label = str(label) if label else None
        cls.value = float(value) if value else None
        cls.unit = str(unit).split("/")[-1] if unit else None
        return cls

    @classmethod
    def from_graph(cls, graph, id):
        cls = cls._get_value(graph, id)
        return cls


class InputParameter(Property):
    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.InputParameter, label=self.label)
        self._add_value(graph, property)


class OutputParameter(Property):
    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.OutputParameter, label=self.label)
        self._add_value(graph, property)


class CalculatedProperty(Property):
    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.CalculatedProperty, label=self.label)
        self._add_value(graph, property)


class EnergyCutoff(InputParameter):
    pid: str = ASMO.EnergyCutoff.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.EnergyCutoff, label=self.label)
        self._add_value(graph, property)


class KPointMesh(InputParameter):
    pid: str = ASMO.KPointMesh.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.KPointMesh, label=self.label)
        self._add_value(graph, property)


class NumberOfIonicSteps(InputParameter):
    pid: str = ASMO.NumberOfIonicSteps.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.NumberOfIonicSteps, label=self.label)
        self._add_value(graph, property)


class PeriodicBoundaryCondition(InputParameter):
    pid: str = ASMO.PeriodicBoundaryCondition.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(
            name, ASMO.PeriodicBoundaryCondition, label=self.label
        )
        self._add_value(graph, property)


class TimeStep(InputParameter):
    pid: str = ASMO.TimeStep.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.TimeStep, label=self.label)
        self._add_value(graph, property)


class VolumeRange(InputParameter):
    pid: str = ASMO.VolumeRange.uri

    def to_graph(self, graph):
        # this knows just the serialisation of itself.
        name = self._create_name()
        self.id = name
        property = graph.create_node(name, ASMO.VolumeRange, label=self.label)
        self._add_value(graph, property)
