from typing import List, Optional, Union
from atomrdf import graph
from pydantic import BaseModel, Field, field_validator
import uuid
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
from atomrdf.datamodels.structure import AtomicScaleSample
from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.utils import get_sample_object
from ase import Atoms
from atomrdf.build.bulk import _generate_atomic_sample_data
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

    @field_validator("initial_sample", mode="before")
    @classmethod
    def _validate_initial_sample(cls, sample):
        if isinstance(sample, Atoms):
            if "id" in sample.info:
                return sample.info["id"]
            else:
                #we have serialise the structure and then return the id
                data = _generate_atomic_sample_data(sample)
                sample_obj = AtomicScaleSample(**data)
                return sample_obj
        elif isinstance(sample, AtomicScaleSample):
            return sample
        elif isinstance(sample, str):
            return sample

    @field_validator("final_sample", mode="before")
    @classmethod
    def _validate_final_sample(cls, sample):
        if isinstance(sample, Atoms):
            if "id" in sample.info:
                return sample.info["id"]
            else:
                #we have serialise the structure and then return the id
                data = _generate_atomic_sample_data(sample)
                sample_obj = AtomicScaleSample(**data)
                return sample_obj
        elif isinstance(sample, AtomicScaleSample):
            return sample
        elif isinstance(sample, str):
            return sample
