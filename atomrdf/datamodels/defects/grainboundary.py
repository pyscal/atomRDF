from typing import List, Optional, Union
import numpy as np
from atomrdf import graph
from pydantic import Field
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty, BaseModel
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
from atomrdf.utils import get_material


class GrainBoundary(TemplateMixin, BaseModel):
    sigma: Optional[int] = None
    plane: Optional[List[float]] = None
    misorientation_angle: Optional[float] = None
    rotation_axis: Optional[List[float]] = None

    def _add_gb(self, graph, name, plane_defect):
        graph.add(
            (
                plane_defect,
                PLDO.hasSigmaValue,
                Literal(self.sigma, datatype=XSD.integer),
            )
        )
        graph.add(
            (
                plane_defect,
                PLDO.hasGBplane,
                Literal(
                    " ".join(np.array(self.plane).astype(str)), datatype=XSD.string
                ),
            )
        )
        graph.add(
            (
                plane_defect,
                PLDO.hasRotationAxis,
                Literal(
                    " ".join(np.array(self.rotation_axis).astype(str)),
                    datatype=XSD.string,
                ),
            )
        )
        graph.add(
            (
                plane_defect,
                PLDO.hasMisorientationAngle,
                Literal(self.misorientation_angle, datatype=XSD.float),
            )
        )

    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        plane_defect = graph.create_node(f"{name}_GrainBoundary", PLDO.GrainBoundary)
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)

    @classmethod
    def _read_gb(cls, graph, sample, plane_defect):
        sigma = graph.value(plane_defect, PLDO.hasSigmaValue)
        plane = graph.value(plane_defect, PLDO.hasGBplane)
        rotation_axis = graph.value(plane_defect, PLDO.hasRotationAxis)
        misorientation_angle = graph.value(plane_defect, PLDO.hasMisorientationAngle)

        return cls(
            sigma=int(sigma) if sigma is not None else None,
            plane=[float(x) for x in plane.split()] if plane is not None else None,
            rotation_axis=(
                [float(x) for x in rotation_axis.split()]
                if rotation_axis is not None
                else None
            ),
            misorientation_angle=(
                float(misorientation_angle)
                if misorientation_angle is not None
                else None
            ),
        )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            plane_defect = triple[2]
            typev = graph.value(plane_defect, RDF.type)
            if typev is not None and typev.toPython() == PLDO.GrainBoundary.uri:
                return cls._read_gb(graph, sample, plane_defect)
        return None


class TwistGrainBoundary(GrainBoundary):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        plane_defect = graph.create_node(
            f"{name}_TwistGrainBoundary", PLDO.TwistGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            plane_defect = triple[2]
            typev = graph.value(plane_defect, RDF.type)
            if typev is not None and typev.toPython() == PLDO.TwistGrainBoundary.uri:
                return cls._read_gb(graph, sample, plane_defect)
        return None


class TiltGrainBoundary(GrainBoundary):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        plane_defect = graph.create_node(
            f"{name}_TiltGrainBoundary", PLDO.TiltGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            plane_defect = triple[2]
            typev = graph.value(plane_defect, RDF.type)
            if typev is not None and typev.toPython() == PLDO.TiltGrainBoundary.uri:
                return cls._read_gb(graph, sample, plane_defect)
        return None


class SymmetricalTiltGrainBoundary(GrainBoundary):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        print("point 1")
        plane_defect = graph.create_node(
            f"{name}_SymmetricalTiltGrainBoundary", PLDO.SymmetricalTiltGrainBoundary
        )
        print("point 2")
        print(plane_defect)
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            plane_defect = triple[2]
            typev = graph.value(plane_defect, RDF.type)
            if (
                typev is not None
                and typev.toPython() == PLDO.SymmetricalTiltGrainBoundary.uri
            ):
                return cls._read_gb(graph, sample, plane_defect)
        return None


class MixedGrainBoundary(GrainBoundary):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        plane_defect = graph.create_node(
            f"{name}_MixedGrainBoundary", PLDO.MixedGrainBoundary
        )
        graph.add((material, CDCO.hasCrystallographicDefect, plane_defect))
        self._add_gb(graph, name)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            plane_defect = triple[2]
            typev = graph.value(plane_defect, RDF.type)
            if typev is not None and typev.toPython() == PLDO.MixedGrainBoundary.uri:
                return cls._read_gb(graph, sample, plane_defect)
        return None
