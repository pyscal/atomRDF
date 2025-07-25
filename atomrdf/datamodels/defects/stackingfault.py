from typing import List, Optional, Union
import numpy as np
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


class StackingFault(TemplateMixin, BaseModel):
    plane: Optional[DataProperty[List[float]]] = None
    displacement: Optional[DataProperty[List[float]]] = None

    def to_graph(self, graph, name, material):
        plane = " ".join(np.array(self.plane.value).astype(str))
        displ = " ".join(np.array(self.displacement.value).astype(str))
        sf = graph.create_node(f"{name}_StackingFault", PLDO.StackingFault)
        graph.add((material, CDCO.hasCrystallographicDefect, sf))
        graph.add((sf, PLDO.hasSFplane, Literal(plane, datatype=XSD.string)))
        graph.add((sf, PLDO.hasDisplacementVector, Literal(displ, datatype=XSD.string)))
