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


class PointDefect(TemplateMixin, BaseModel):
    concentration: Optional[DataProperty[float]] = None
    number: Optional[DataProperty[int]] = None


class Vacancy(PointDefect):
    def to_graph(self, graph, name, material, sample):
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


class Substitutional(PointDefect):
    def to_graph(self, graph, name, material, sample):
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


class Interstitial(PointDefect):
    def to_graph(self, graph, name, material, sample):
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
