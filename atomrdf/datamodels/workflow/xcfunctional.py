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


class XCFunctional(TemplateMixin, BaseModel):
    pid: str = MDO.ExchangeCorrelationEnergyFunctional.uri

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return MDO.ExchangeCorrelationEnergyFunctional


class GGA(XCFunctional):
    pid: str = MDO.GGA.uri

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return MDO.GGA


class LDA(XCFunctional):
    pid: str = MDO.LDA.uri

    def to_graph(self):
        # as strange as it may seem, this is what a NamedIndividual should do
        return MDO.LDA
