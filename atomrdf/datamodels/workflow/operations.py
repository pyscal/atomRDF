from typing import List, Optional, Union
from atomrdf import graph
from pydantic import BaseModel, Field
import uuid
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty, Activity
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
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        final_sample = graph.value(activity, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)
        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
        )


class SubstituteAtom(Activity):
    def to_graph(self, graph):
        activity_id = f"substituteatom:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.SubstituteAtom)
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        final_sample = graph.value(activity, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)
        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
        )


class AddAtom(Activity):
    def to_graph(self, graph):
        activity_id = f"addatom:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.AddAtom)
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        final_sample = graph.value(activity, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)
        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
        )


class Rotate(Activity):
    rotation_matrix: Optional[DataProperty[List[List[float]]]] = None

    def to_graph(self, graph):
        activity_id = f"rotate:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Rotation)
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

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
        final_sample = graph.value(activity_id, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)

        rotation_matrix = []
        rot_matrix = graph.objects(activity, CMSO.hasVector)
        for rot_vector in rot_matrix:
            x = graph.value(rot_vector, CMSO.hasComponent_x)
            y = graph.value(rot_vector, CMSO.hasComponent_y)
            z = graph.value(rot_vector, CMSO.hasComponent_z)
            rotation_matrix.append([x.toPython(), y.toPython(), z.toPython()])

        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
            rotation_matrix=DataProperty(value=rotation_matrix),
        )


class Translate(Activity):
    translation_vector: Optional[DataProperty[List[float]]] = None

    def to_graph(self, graph):
        activity_id = f"translate:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Translation)
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

        translation_vector = graph.create_node(
            f"{activity_id}_TranslationVector", CMSO.Vector
        )
        graph.add((activity, CMSO.hasVector, translation_vector))
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_x,
                Literal(self.translation_vector.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_y,
                Literal(self.translation_vector.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                translation_vector,
                CMSO.hasComponent_z,
                Literal(self.translation_vector.value[2], datatype=XSD.float),
            )
        )

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        final_sample = graph.value(activity_id, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)

        translation_vector_lst = []
        translation_vector = graph.objects(activity, CMSO.hasVector)
        x = graph.value(translation_vector, CMSO.hasComponent_x)
        y = graph.value(translation_vector, CMSO.hasComponent_y)
        z = graph.value(translation_vector, CMSO.hasComponent_z)
        translation_vector_lst.append([x.toPython(), y.toPython(), z.toPython()])
        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
            translation_vector=DataProperty(value=translation_vector_lst),
        )


class Shear(Activity):
    shear_vector: Optional[DataProperty[List[float]]] = None
    normal_vector: Optional[DataProperty[List[float]]] = None

    def to_graph(self, graph):
        activity_id = f"shear:{str(uuid.uuid4())}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Shear)
        graph.add((self.final_sample, PROV.wasDerivedFrom, self.initial_sample))
        graph.add((self.final_sample, PROV.wasGeneratedBy, activity))

        shear_vector = graph.create_node(f"{activity_id}_ShearVector", CMSO.Vector)
        graph.add((activity, CMSO.hasVector, shear_vector))
        graph.add(
            (
                shear_vector,
                CMSO.hasComponent_x,
                Literal(self.shear_vector.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                shear_vector,
                CMSO.hasComponent_y,
                Literal(self.shear_vector.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                shear_vector,
                CMSO.hasComponent_z,
                Literal(self.shear_vector.value[2], datatype=XSD.float),
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
                    Literal(self.normal_vector.value[0], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    normal_vector,
                    CMSO.hasComponent_y,
                    Literal(self.normal_vector.value[1], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    normal_vector,
                    CMSO.hasComponent_z,
                    Literal(self.normal_vector.value[2], datatype=XSD.float),
                )
            )

    @classmethod
    def from_graph(cls, graph, activity_id):
        activity = get_sample_object(activity_id)
        final_sample = graph.value(activity_id, PROV.wasGeneratedBy)
        initial_sample = graph.value(final_sample, PROV.wasDerivedFrom)

        shear_vector_lst = []
        shear_vector = graph.objects(activity, CMSO.hasVector)
        x = graph.value(shear_vector, CMSO.hasComponent_x)
        y = graph.value(shear_vector, CMSO.hasComponent_y)
        z = graph.value(shear_vector, CMSO.hasComponent_z)
        shear_vector_lst.append([x.toPython(), y.toPython(), z.toPython()])

        normal_vector_lst = []
        plane = graph.objects(activity, CMSO.hasPlane)
        if plane:
            normal_vector = graph.objects(plane, CMSO.hasNormalVector)
            x = graph.value(normal_vector, CMSO.hasComponent_x)
            y = graph.value(normal_vector, CMSO.hasComponent_y)
            z = graph.value(normal_vector, CMSO.hasComponent_z)
            normal_vector_lst.append([x.toPython(), y.toPython(), z.toPython()])
        if len(normal_vector_lst) == 0:
            normal_vector_lst = None

        return cls(
            id=activity_id,
            initial_sample=initial_sample,
            final_sample=final_sample,
            shear_vector=DataProperty(value=shear_vector_lst),
            normal_vector=DataProperty(value=normal_vector_lst),
        )
