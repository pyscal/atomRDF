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


class PointDefect(TemplateMixin, BaseModel):
    concentration: Optional[DataProperty[float]] = None
    number: Optional[DataProperty[int]] = None


class Vacancy(PointDefect):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        sample = URIRef(sample_id)
        vacancy = graph.create_node(f"{name}_Vacancy", PODO.Vacancy)
        graph.add((material, CDCO.hasCrystallographicDefect, vacancy))
        graph.add(
            (
                sample,
                PODO.hasVacancyConcentration,
                Literal(self.concentration.value, datatype=XSD.float),
            )
        )
        if self.number is not None:
            graph.add(
                (
                    sample,
                    PODO.hasNumberOfVacancies,
                    Literal(self.number.value, datatype=XSD.integer),
                )
            )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            vacancy = triple[2]
            typev = graph.value(vacancy, RDF.type)
            if typev is not None and typev.toPython() == PODO.Vacancy.uri:
                concentration = graph.value(sample, PODO.hasVacancyConcentration)
                number = graph.value(sample, PODO.hasNumberOfVacancies)
                return cls(
                    concentration=DataProperty(value=concentration),
                    number=DataProperty(value=number),
                )


class Substitutional(PointDefect):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        sample = URIRef(sample_id)
        defect = graph.create_node(
            f"{name}_SubstitutionalImpurity", PODO.SubstitutionalImpurity
        )
        graph.add((material, CDCO.hasCrystallographicDefect, defect))
        graph.add(
            (
                sample,
                PODO.hasImpurityConcentration,
                Literal(self.concentration.value, datatype=XSD.float),
            )
        )
        if self.number is not None:
            graph.add(
                (
                    sample,
                    PODO.hasNumberOfImpurityAtoms,
                    Literal(self.number.value, datatype=XSD.integer),
                )
            )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            defect = triple[2]
            typev = graph.value(defect, RDF.type)
            if (
                typev is not None
                and typev.toPython() == PODO.SubstitutionalImpurity.uri
            ):
                concentration = graph.value(sample, PODO.hasImpurityConcentration)
                number = graph.value(sample, PODO.hasNumberOfImpurityAtoms)
                return cls(
                    concentration=DataProperty(value=concentration),
                    number=DataProperty(value=number),
                )


class Interstitial(PointDefect):
    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        sample = URIRef(sample_id)
        defect = graph.create_node(
            f"{name}_InterstitialImpurity", PODO.InterstitialImpurity
        )
        graph.add((material, CDCO.hasCrystallographicDefect, defect))
        graph.add(
            (
                sample,
                PODO.hasImpurityConcentration,
                Literal(self.concentration.value, datatype=XSD.float),
            )
        )
        if self.number is not None:
            graph.add(
                (
                    sample,
                    PODO.hasNumberOfImpurityAtoms,
                    Literal(self.number.value, datatype=XSD.integer),
                )
            )

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            defect = triple[2]
            typev = graph.value(defect, RDF.type)
            if typev is not None and typev.toPython() == PODO.InterstitialImpurity.uri:
                concentration = graph.value(sample, PODO.hasImpurityConcentration)
                number = graph.value(sample, PODO.hasNumberOfImpurityAtoms)
                return cls(
                    concentration=DataProperty(value=concentration),
                    number=DataProperty(value=number),
                )
