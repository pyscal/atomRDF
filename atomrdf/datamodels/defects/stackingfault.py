from typing import List, Optional, Union
import numpy as np
from atomrdf import graph
from pydantic import Field
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty, BaseModel
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


class StackingFault(TemplateMixin, BaseModel):
    plane: Optional[List[float]] = None
    displacement: Optional[List[float]] = None

    def to_graph(self, graph, sample_id):
        name = sample_id
        material = get_material(graph, sample_id)
        plane = " ".join(np.array(self.plane).astype(str))
        displ = " ".join(np.array(self.displacement).astype(str))
        sf = graph.create_node(f"{name}_StackingFault", PLDO.StackingFault)
        graph.add((material, CDCO.hasCrystallographicDefect, sf))
        graph.add((sf, PLDO.hasSFplane, Literal(plane, datatype=XSD.string)))
        graph.add((sf, PLDO.hasDisplacementVector, Literal(displ, datatype=XSD.string)))

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        for triple in graph.triples((material, CDCO.hasCrystallographicDefect, None)):
            sf = triple[2]
            typev = graph.value(sf, RDF.type)
            if typev is not None and typev.toPython() == PLDO.StackingFault.uri:
                plane = graph.value(sf, PLDO.hasSFplane)
                displacement = graph.value(sf, PLDO.hasDisplacementVector)
                return cls(
                    plane=(
                        [float(x) for x in plane.toPython().split()]
                        if plane is not None
                        else None
                    ),
                    displacement=(
                        [float(x) for x in displacement.toPython().split()]
                        if displacement is not None
                        else None
                    ),
                )
        return None
