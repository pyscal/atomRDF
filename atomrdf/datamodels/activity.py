from typing import List, Optional, Union, TYPE_CHECKING
from pydantic import BaseModel, Field, field_validator
import uuid
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
from atomrdf.datamodels.structure import AtomicScaleSample

if TYPE_CHECKING:
    from atomrdf.graph import KnowledgeGraph
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
    model_config = {"arbitrary_types_allowed": True}

    pid: Optional[str] = Field(default=None, description="PID of the activity")

    # Graph reference for loading samples from strings
    graph: Optional["KnowledgeGraph"] = Field(
        default=None, exclude=True, description="Knowledge graph reference"
    )

    input_sample: Optional[
        Union[str, AtomicScaleSample, List[Union[str, AtomicScaleSample]]]
    ] = Field(
        default=None,
        description="Input sample(s) - can be AtomicScaleSample object, ID string, or list of either",
    )
    output_sample: Optional[
        Union[str, AtomicScaleSample, List[Union[str, AtomicScaleSample]]]
    ] = Field(
        default=None,
        description="Output sample(s) - can be AtomicScaleSample object, ID string, or list of either",
    )

    @field_validator("input_sample", mode="before")
    @classmethod
    def _validate_input_sample(cls, sample, info):
        # Handle list of samples
        if isinstance(sample, list):
            return [cls._validate_single_sample(s, info) for s in sample]
        return cls._validate_single_sample(sample, info)

    @staticmethod
    def _validate_single_sample(sample, info):
        if isinstance(sample, Atoms):
            if "id" in sample.info:
                return sample.info["id"]
            else:
                # We have to serialize the structure and then return the object
                data = _generate_atomic_sample_data(sample)
                sample_obj = AtomicScaleSample(**data)
                return sample_obj
        elif isinstance(sample, AtomicScaleSample):
            return sample
        elif isinstance(sample, str):
            # If string and graph is available, load from graph
            graph_obj = info.data.get("graph") if info.data else None
            if graph_obj is not None:
                try:
                    return AtomicScaleSample.from_graph(graph_obj, sample)
                except:
                    # If loading fails, keep as string
                    return sample
            return sample
        return sample

    @field_validator("output_sample", mode="before")
    @classmethod
    def _validate_output_sample(cls, sample, info):
        # Handle list of samples
        if isinstance(sample, list):
            return [cls._validate_single_sample(s, info) for s in sample]
        return cls._validate_single_sample(sample, info)
