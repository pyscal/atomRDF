from typing import List, Optional, Union
from atomrdf import graph
from pydantic import BaseModel, Field
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
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


class GrainBoundary(TemplateMixin, BaseModel):
    sigma: Optional[DataProperty[int]] = None
    plane: Optional[DataProperty[List[float]]] = None
    misorientation_angle: Optional[DataProperty[float]] = None
    rotation_axis: Optional[DataProperty[List[float]]] = None

    def _add_gb(self, graph, name, plane_defect):
        graph.add(
            (
                plane_defect,
                PLDO.hasSigmaValue,
                Literal(self.sigma.value, datatype=XSD.integer),
            )
        )
        graph.add(
            (
                plane_defect,
                PLDO.hasGBplane,
                Literal(" ".join(self.plane.value.astype(str)), datatype=XSD.string),
            )
        )
        graph.add(
            (
                plane_defect,
                PLDO.hasRotationAxis,
                Literal(
                    " ".join(self.rotation_axis.value.astype(str)), datatype=XSD.string
                ),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasMisorientationAngle,
                Literal(self.misorientation_angle.value, datatype=XSD.float),
            )
        )

    def to_graph(self, graph, name, material):
        plane_defect = graph.create_node(f"{name}_GrainBoundary", PLDO.GrainBoundary)
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)


class TwistGrainBoundary(GrainBoundary):
    def to_graph(self, graph, name, material):
        plane_defect = graph.create_node(
            f"{name}_TwistGrainBoundary", PLDO.TwistGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)


class TiltGrainBoundary(GrainBoundary):
    def to_graph(self, graph, name, material):
        plane_defect = graph.create_node(
            f"{name}_TiltGrainBoundary", PLDO.TiltGrainBoundary
        )
        self._add_gb(graph, name)


class SymmetricalTiltGrainBoundary(GrainBoundary):
    def to_graph(self, graph, name, material):
        plane_defect = graph.create_node(
            f"{name}_SymmetricalTiltGrainBoundary", PLDO.SymmetricalTiltGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)


class MixedGrainBoundary(GrainBoundary):
    def to_graph(self, graph, name, material):
        plane_defect = graph.create_node(
            f"{name}_MixedGrainBoundary", PLDO.MixedGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)
