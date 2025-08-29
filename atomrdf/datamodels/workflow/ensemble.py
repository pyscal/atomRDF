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


class ThermodynamicEnsemble(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the ensemble")


class CanonicalEnsemble(ThermodynamicEnsemble):
    pid: str = ASMO.CanonicalEnsemble.uri
    basename: str = "CanonicalEnsemble"


class MicrocanonicalEnsemble(ThermodynamicEnsemble):
    pid: str = ASMO.MicrocanonicalEnsemble.uri
    basename: str = "MicrocanonicalEnsemble"


class IsothermalIsobaricEnsemble(ThermodynamicEnsemble):
    pid: str = ASMO.IsothermalIsobaricEnsemble.uri
    basename: str = "IsothermalIsobaricEnsemble"


class IsoenthalpicIsobaricEnsemble(ThermodynamicEnsemble):
    pid: str = ASMO.IsoenthalpicIsobaricEnsemble.uri
    basename: str = "IsoenthalpicIsobaricEnsemble"


class GrandCanonicalEnsemble(ThermodynamicEnsemble):
    pid: str = ASMO.GrandCanonicalEnsemble.uri
    basename: str = "GrandCanonicalEnsemble"
