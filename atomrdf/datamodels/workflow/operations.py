from typing import List, Optional, Union
from atomrdf import graph
from pydantic import BaseModel, Field
import uuid
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
from atomrdf.datamodels.activity import Activity
from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.utils import get_sample_object
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


class DeleteAtom(Activity):
    def to_graph(self, graph):
        activity_id = f"deleteatom:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.DeleteAtom)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None
        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
        )


class SubstituteAtom(Activity):
    def to_graph(self, graph):
        activity_id = f"substituteatom:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.SubstituteAtom)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None
        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
        )


class AddAtom(Activity):
    def to_graph(self, graph):
        activity_id = f"addatom:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.AddAtom)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None
        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
        )


class Rotate(Activity):
    rotation_matrix: Optional[List[List[float]]] = None

    def to_graph(self, graph):
        activity_id = f"rotate:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Rotation)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

        rot_vector_01 = graph.create_node(
            f"{activity_id}_RotationVector_1", CMSO.Vector
        )
        graph.add((activity, CMSO.hasVector, rot_vector_01))
        graph.add(
            (
                rot_vector_01,
                CMSO.hasComponent_x,
                Literal(self.rotation_matrix[0][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_01,
                CMSO.hasComponent_y,
                Literal(self.rotation_matrix[0][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_01,
                CMSO.hasComponent_z,
                Literal(self.rotation_matrix[0][2], datatype=XSD.float),
            )
        )

        rot_vector_02 = graph.create_node(
            f"{activity_id}_RotationVector_2", CMSO.Vector
        )
        graph.add((activity, CMSO.hasVector, rot_vector_02))
        graph.add(
            (
                rot_vector_02,
                CMSO.hasComponent_x,
                Literal(self.rotation_matrix[1][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_02,
                CMSO.hasComponent_y,
                Literal(self.rotation_matrix[1][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_02,
                CMSO.hasComponent_z,
                Literal(self.rotation_matrix[1][2], datatype=XSD.float),
            )
        )

        rot_vector_03 = graph.create_node(
            f"{activity_id}_RotationVector_3", CMSO.Vector
        )
        graph.add((activity, CMSO.hasVector, rot_vector_03))
        graph.add(
            (
                rot_vector_03,
                CMSO.hasComponent_x,
                Literal(self.rotation_matrix[2][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_03,
                CMSO.hasComponent_y,
                Literal(self.rotation_matrix[2][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                rot_vector_03,
                CMSO.hasComponent_z,
                Literal(self.rotation_matrix[2][2], datatype=XSD.float),
            )
        )

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None

        rotation_matrix = []
        for rot_vector in graph.objects(activity, CMSO.hasVector):
            x = graph.value(rot_vector, CMSO.hasComponent_x)
            y = graph.value(rot_vector, CMSO.hasComponent_y)
            z = graph.value(rot_vector, CMSO.hasComponent_z)
            rotation_matrix.append([x.toPython(), y.toPython(), z.toPython()])

        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
            rotation_matrix=rotation_matrix,
        )


class Translate(Activity):
    translation_vector: Optional[List[float]] = None

    def to_graph(self, graph):
        activity_id = f"translate:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Translation)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

        translation_vector = graph.create_node(
            f"{activity_id}_TranslationVector", CMSO.Vector
        )
        graph.add((activity, CMSO.hasVector, translation_vector))
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_x,
                Literal(self.translation_vector[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_y,
                Literal(self.translation_vector[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_z,
                Literal(self.translation_vector[2], datatype=XSD.float),
            )
        )

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None

        translation_vector = None
        for vec in graph.objects(activity, CMSO.hasVector):
            x = graph.value(vec, CMSO.hasComponent_x)
            y = graph.value(vec, CMSO.hasComponent_y)
            z = graph.value(vec, CMSO.hasComponent_z)
            translation_vector = [x.toPython(), y.toPython(), z.toPython()]
            break
        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
            translation_vector=translation_vector,
        )


class Shear(Activity):
    shear_vector: Optional[List[float]] = None
    normal_vector: Optional[List[float]] = None
    distance: Optional[float] = None

    def to_graph(self, graph):
        activity_id = f"shear:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Shear)
        graph.add(
            (URIRef(self.output_sample), PROV.wasDerivedFrom, URIRef(self.input_sample))
        )
        graph.add((URIRef(self.output_sample), PROV.wasGeneratedBy, activity))

        if self.shear_vector is not None:
            shear_vector = graph.create_node(f"{activity_id}_ShearVector", CMSO.Vector)
            graph.add((activity, CMSO.hasVector, shear_vector))
            graph.add(
                (
                    shear_vector,
                    CMSO.hasComponent_x,
                    Literal(self.shear_vector[0], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    shear_vector,
                    CMSO.hasComponent_y,
                    Literal(self.shear_vector[1], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    shear_vector,
                    CMSO.hasComponent_z,
                    Literal(self.shear_vector[2], datatype=XSD.float),
                )
            )

        if self.normal_vector:
            plane = graph.create_node(f"{activity_id}_Plane", CMSO.Plane)
            normal_vector = graph.create_node(
                f"{activity_id}_NormalVector", CMSO.NormalVector
            )
            graph.add((activity, CMSO.hasPlane, plane))
            graph.add((plane, CMSO.hasNormalVector, normal_vector))
            graph.add(
                (
                    normal_vector,
                    CMSO.hasComponent_x,
                    Literal(self.normal_vector[0], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    normal_vector,
                    CMSO.hasComponent_y,
                    Literal(self.normal_vector[1], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    normal_vector,
                    CMSO.hasComponent_z,
                    Literal(self.normal_vector[2], datatype=XSD.float),
                )
            )
            if self.distance:
                graph.add(
                    (
                        plane,
                        CMSO.hasDistanceFromOrigin,
                        Literal(self.distance, datatype=XSD.float),
                    )
                )

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        output_sample = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        input_sample = graph.value(output_sample, PROV.wasDerivedFrom) if output_sample else None

        shear_vector = None
        for vec in graph.objects(activity, CMSO.hasVector):
            x = graph.value(vec, CMSO.hasComponent_x)
            y = graph.value(vec, CMSO.hasComponent_y)
            z = graph.value(vec, CMSO.hasComponent_z)
            shear_vector = [x.toPython(), y.toPython(), z.toPython()]
            break

        normal_vector = None
        distance = None
        plane = graph.value(activity, CMSO.hasPlane)
        if plane is not None:
            nv = graph.value(plane, CMSO.hasNormalVector)
            if nv is not None:
                x = graph.value(nv, CMSO.hasComponent_x)
                y = graph.value(nv, CMSO.hasComponent_y)
                z = graph.value(nv, CMSO.hasComponent_z)
                normal_vector = [x.toPython(), y.toPython(), z.toPython()]
            d = graph.value(plane, CMSO.hasDistanceFromOrigin)
            if d is not None:
                distance = d.toPython()

        return cls(
            id=str(activity_id),
            input_sample=str(input_sample) if input_sample else None,
            output_sample=str(output_sample) if output_sample else None,
            shear_vector=shear_vector,
            normal_vector=normal_vector,
            distance=distance,
        )
