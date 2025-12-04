from typing import List, Optional, Union
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


class SlipPlane(TemplateMixin, BaseModel):
    normal: Optional[List[float]] = None

    def to_graph(self, graph, sample_id, slip_system):
        slip_plane = graph.create_node(
            f"{sample_id}_DislocationSlipPlane", LDO.SlipPlane
        )
        normal_vector = graph.create_node(
            f"{sample_id}_DislocationNormalVector", LDO.NormalVector
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_x,
                Literal(self.normal[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_y,
                Literal(self.normal[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                normal_vector,
                CMSO.hasComponent_z,
                Literal(self.normal[2], datatype=XSD.float),
            )
        )
        graph.add((slip_plane, LDO.hasNormalVector, normal_vector))
        graph.add((slip_plane, LDO.belongsToSystem, slip_system))

    @classmethod
    def from_graph(cls, graph, sample, slip_system):
        for triple in graph.triples((None, LDO.belongsToSystem, slip_system)):
            typev = graph.value(triple[0], RDF.type)
            if typev is not None:
                if typev.toPython() == LDO.SlipPlane.uri:
                    slip_plane = triple[0]
                    break

        normal_vector = graph.value(slip_plane, LDO.hasNormalVector)
        normal_x = graph.value(normal_vector, CMSO.hasComponent_x)
        normal_y = graph.value(normal_vector, CMSO.hasComponent_y)
        normal_z = graph.value(normal_vector, CMSO.hasComponent_z)

        return cls(normal=[float(normal_x), float(normal_y), float(normal_z)])


class SlipSystem(TemplateMixin, BaseModel):
    slip_direction: Optional[List[float]] = None
    slip_plane: Optional[SlipPlane] = None

    def to_graph(self, graph, sample_id, line_defect):
        slip_system = graph.create_node(
            f"{sample_id}_DislocationSlipSystem", LDO.SlipSystem
        )
        slip_direction = graph.create_node(
            f"{sample_id}_DislocationSlipDirection", LDO.SlipDirection
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_x,
                Literal(self.slip_direction[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_y,
                Literal(self.slip_direction[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                slip_direction,
                CMSO.hasComponent_z,
                Literal(self.slip_direction[2], datatype=XSD.float),
            )
        )
        graph.add((slip_direction, LDO.belongsToSystem, slip_system))
        graph.add((line_defect, LDO.movesOn, slip_system))
        self.slip_plane.to_graph(graph, sample_id, slip_system)

    @classmethod
    def from_graph(cls, graph, sample, line_defect):
        slip_system = graph.value(line_defect, LDO.movesOn)
        for triple in graph.triples((None, LDO.belongsToSystem, slip_system)):
            typev = graph.value(triple[0], RDF.type)
            if typev is not None:
                if typev.toPython() == LDO.SlipDirection.uri:
                    slip_direction = triple[0]
                    break

        direction_x = graph.value(slip_direction, CMSO.hasComponent_x)
        direction_y = graph.value(slip_direction, CMSO.hasComponent_y)
        direction_z = graph.value(slip_direction, CMSO.hasComponent_z)

        slip_plane_instance = SlipPlane.from_graph(graph, sample, slip_system)

        return cls(
            slip_direction=[float(direction_x), float(direction_y), float(direction_z)],
            slip_plane=slip_plane_instance,
        )


class Dislocation(TemplateMixin, BaseModel):
    line_direction: Optional[List[float]] = None
    burgers_vector: Optional[List[float]] = None
    slip_system: Optional[SlipSystem] = None

    def _add_dislocation(self, graph, sample_id, line_defect):
        line_direction = graph.create_node(
            f"{sample_id}_DislocationLineDirection", LDO.LineDirection
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_x,
                Literal(self.line_direction[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_y,
                Literal(self.line_direction[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                line_direction,
                CMSO.hasComponent_z,
                Literal(self.line_direction[2], datatype=XSD.float),
            )
        )
        graph.add((line_defect, LDO.hasLineDirection, line_direction))

        burgers_vector = graph.create_node(
            f"{sample_id}_DislocationBurgersVector", LDO.BurgersVector
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_x,
                Literal(self.burgers_vector[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_y,
                Literal(self.burgers_vector[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                burgers_vector,
                CMSO.hasComponent_z,
                Literal(self.burgers_vector[2], datatype=XSD.float),
            )
        )
        graph.add((line_defect, LDO.hasBurgersVector, burgers_vector))
        self.slip_system.to_graph(graph, sample_id, line_defect)

    def to_graph(self, graph, sample_id):
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{sample_id}_Dislocation", LDO.Dislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, sample_id, line_defect)

    @classmethod
    def _read_dislocation(cls, graph, sample, line_defect):
        line_direction = graph.value(line_defect, LDO.hasLineDirection)
        burgers_vector = graph.value(line_defect, LDO.hasBurgersVector)

        direction_x = graph.value(line_direction, CMSO.hasComponent_x)
        direction_y = graph.value(line_direction, CMSO.hasComponent_y)
        direction_z = graph.value(line_direction, CMSO.hasComponent_z)

        burgers_x = graph.value(burgers_vector, CMSO.hasComponent_x)
        burgers_y = graph.value(burgers_vector, CMSO.hasComponent_y)
        burgers_z = graph.value(burgers_vector, CMSO.hasComponent_z)

        slip_system_instance = SlipSystem.from_graph(graph, sample, line_defect)

        return cls(
            line_direction=[float(direction_x), float(direction_y), float(direction_z)],
            burgers_vector=[float(burgers_x), float(burgers_y), float(burgers_z)],
            slip_system=slip_system_instance,
        )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            line_defect = triple[2]
            typev = graph.value(line_defect, RDF.type)
            if typev is not None and typev.toPython() == LDO.Dislocation.uri:
                return cls._read_dislocation(graph, sample, line_defect)


class ScrewDislocation(Dislocation):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.ScrewDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            line_defect = triple[2]
            typev = graph.value(line_defect, RDF.type)
            if typev is not None and typev.toPython() == LDO.ScrewDislocation.uri:
                return cls._read_dislocation(graph, sample, line_defect)


class EdgeDislocation(Dislocation):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.EdgeDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            line_defect = triple[2]
            typev = graph.value(line_defect, RDF.type)
            if typev is not None and typev.toPython() == LDO.EdgeDislocation.uri:
                return cls._read_dislocation(graph, sample, line_defect)


class MixedDislocation(Dislocation):
    character_angle: Optional[float] = None

    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        line_defect = graph.create_node(f"{name}_Dislocation", LDO.MixedDislocation)
        graph.add((material, CDCO.hasCrystallographicDefect, line_defect))
        self._add_dislocation(graph, name, line_defect)
        graph.add(
            (
                line_defect,
                LDO.hasCharacterAngle,
                Literal(self.character_angle, datatype=XSD.float),
            )
        )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            line_defect = triple[2]
            typev = graph.value(line_defect, RDF.type)
            if typev is not None and typev.toPython() == LDO.MixedDislocation.uri:
                dislocation = cls._read_dislocation(graph, sample, line_defect)
                character_angle = graph.value(line_defect, LDO.hasCharacterAngle)
                dislocation.character_angle = (
                    float(character_angle) if character_angle is not None else None
                )
                return dislocation
