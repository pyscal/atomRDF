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


class XCFunctional(TemplateMixin, BaseModel):
    pid: str = MDO.ExchangeCorrelationEnergyFunctional.toPython()

    def to_graph(self, graph, main_id):
        # as strange as it may seem, this is what a NamedIndividual should do
        main_id = f"{main_id}_xcfunctional"
        potential = graph.create_node(main_id, MDO.ExchangeCorrelationEnergyFunctional)
        return potential

    @classmethod
    def from_graph(cls, graph, potential_id):
        label = graph.value(potential_id, RDFS.label)
        pot_type = graph.value(potential_id, RDF.type)
        if pot_type is not None:
            pot_type = pot_type.toPython().split("/")[-1]
            return cls(label=label, potential_type=pot_type)
        return cls(label=label)


class GGA(XCFunctional):
    pid: str = MDO.GeneralizedGradientApproximation.toPython()

    def to_graph(self, graph, main_id):
        # as strange as it may seem, this is what a NamedIndividual should do
        main_id = f"{main_id}_gga"
        potential = graph.create_node(main_id, MDO.GeneralizedGradientApproximation)
        return potential


class LDA(XCFunctional):
    pid: str = MDO.LocalDensityApproximation.toPython()

    def to_graph(self, graph, main_id):
        # as strange as it may seem, this is what a NamedIndividual should do
        main_id = f"{main_id}_lda"
        potential = graph.create_node(main_id, MDO.LocalDensityApproximation)
        return potential


class HybridFunctional(XCFunctional):
    pid: str = MDO.HybridFunctional.toPython()

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_hybrid_functional"
        potential = graph.create_node(main_id, MDO.HybridFunctional)
        return potential


class HybridGeneralizedGradientApproximation(HybridFunctional):
    pid: str = MDO.HybridGeneralizedGradientApproximation.toPython()

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_hybrid_gga"
        potential = graph.create_node(
            main_id, MDO.HybridGeneralizedGradientApproximation
        )
        return potential


class HybridMetaGeneralizedGradientApproximation(HybridFunctional):
    pid: str = MDO.HybridMetaGeneralizedGradientApproximation.toPython()

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_hybrid_mgga"
        potential = graph.create_node(
            main_id, MDO.HybridMetaGeneralizedGradientApproximation
        )
        return potential


class MetaGeneralizedGradientApproximation(XCFunctional):
    pid: str = MDO.MetaGeneralizedGradientApproximation.toPython()

    def to_graph(self, graph, main_id):
        main_id = f"{main_id}_mgga"
        potential = graph.create_node(main_id, MDO.MetaGeneralizedGradientApproximation)
        return potential


xc_map = {
    "LDA": LDA,
    "GGA": GGA,
    "PBE": GGA,
    "LocalDensityApproximation": LDA,
    "GeneralizedGradientApproximation": GGA,
    "PerdewBurkeErnzerhof": GGA,
    "HybridFunctional": HybridFunctional,
    "HybridGeneralizedGradientApproximation": HybridGeneralizedGradientApproximation,
    "HybridMetaGeneralizedGradientApproximation": HybridMetaGeneralizedGradientApproximation,
    "MetaGeneralizedGradientApproximation": MetaGeneralizedGradientApproximation,
}
