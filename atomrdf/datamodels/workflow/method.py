from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field
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
)


class Method(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the method")
    label: Optional[str] = Field(default=None, description="Label of the method")


class MolecularStatics(Method):
    pid: str = ASMO.MolecularStatics.uri


class MolecularDynamics(Method):
    pid: str = ASMO.MolecularDynamics.uri


class DensityFunctionalTheory(Method):
    pid: str = ASMO.DensityFunctionalTheory.uri


class EquationOfStateFit(Method):
    pid: str = ASMO.EquationOfStateFit.uri


class QuasiHarmonicApproximation(Method):
    pid: str = ASMO.QuasiHarmonicApproximation.uri


class ThermodynamicIntegration(Method):
    pid: str = ASMO.ThermodynamicIntegration.uri
