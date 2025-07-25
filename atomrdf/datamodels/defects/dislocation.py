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
from atomrdf.utils import get_material


class SlipPlane(TemplateMixin, BaseModel):
    normal: Optional[DataProperty[List[float]]] = None

    def to_graph(self, graph, name, slip_system):
        slip_plane = graph.create_node(f"{name}_DislocationSlipPlane", LDO.SlipPlane)
        normal_vector = graph.create_node(
            f"{name}_DislocationNormalVector", LDO.NormalVector
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_x,
                Literal(self.normal.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_y,
                Literal(self.normal.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_z,
                Literal(self.normal.value[2], datatype=XSD.float),
            )
        )
        graph.add((slip_plane, LDO.hasNormalVector, normal_vector))
        graph.add((slip_plane, LDO.belongsToSystem, slip_system))


class SlipSystem(TemplateMixin, BaseModel):
    slip_direction: Optional[DataProperty[List[float]]] = None
    slip_plane: Optional[SlipPlane] = None

    def to_graph(self, graph, name, line_defect):
        slip_system = graph.create_node(f"{name}_DislocationSlipSystem", LDO.SlipSystem)
        slip_direction = graph.create_node(
            f"{name}_DislocationSlipDirection", LDO.SlipDirection
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_x,
                Literal(self.slip_direction.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_y,
                Literal(self.slip_direction.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_z,
                Literal(self.slip_direction.value[2], datatype=XSD.float),
            )
        )
        graph.add((slip_direction, LDO.belongsToSystem, slip_system))
        graph.add((line_defect, LDO.movesOn, slip_system))
        self.slip_plane.to_graph(graph, name, slip_system)


class Dislocation(TemplateMixin, BaseModel):
    line_direction: Optional[DataProperty[List[float]]] = None
    burgers_vector: Optional[DataProperty[List[float]]] = None
    slip_system: Optional[SlipSystem] = None

    def _add_dislocation(self, graph, name, line_defect):
        line_direction = graph.create_node(
            f"{name}_DislocationLineDirection", LDO.LineDirection
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_x,
                Literal(self.line_direction.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_y,
                Literal(self.line_direction.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_z,
                Literal(self.line_direction.value[2], datatype=XSD.float),
            )
        )
        graph.add((line_defect, LDO.hasLineDirection, line_direction))

        burgers_vector = graph.create_node(
            f"{name}_DislocationBurgersVector", LDO.BurgersVector
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_x,
                Literal(self.burgers_vector.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_y,
                Literal(self.burgers_vector.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_z,
                Literal(self.burgers_vector.value[2], datatype=XSD.float),
            )
        )
        graph.add((line_defect, LDO.hasBurgersVector, burgers_vector))
        self.slip_system.to_graph(graph, name, line_defect)

    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.Dislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)


class ScrewDislocation(Dislocation):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.ScrewDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)


class EdgeDislocation(Dislocation):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.EdgeDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)


class MixedDislocation(Dislocation):
    character_angle: Optional[DataProperty[float]] = None

    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.MixedDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)
        self.graph.add(
            (
                line_defect,
                LDO.hasCharacterAngle,
                Literal(self.character_angle.value, datatype=XSD.float),
            )
        )
