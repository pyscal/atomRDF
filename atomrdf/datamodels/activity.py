from typing import List, Optional, Union
from atomrdf import graph
from pydantic import BaseModel, Field
import uuid
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
from atomrdf.datamodels.structure import AtomicScaleSample
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


class Activity(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the activity")

    initial_sample: Optional[Union[str, AtomicScaleSample]] = Field(
        default=None, description="ID of the initial sample in the graph"
    )
    final_sample: Optional[Union[str, AtomicScaleSample]] = Field(
        default=None, description="ID of the final sample in the graph"
    )