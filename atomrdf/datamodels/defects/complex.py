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
    UNSAFECDCO,
)
from atomrdf.utils import get_material


class DefectComplex(TemplateMixin, BaseModel):
    ids: Optional[List[str]] = Field(
        default=None,
        description="List of defect type names (e.g., 'stacking_fault', 'substitutional')",
    )
    relative_distance: Optional[str] = Field(
        default=None, description="Relative distance between defects"
    )

    def to_graph(self, graph, sample_id, defect_uris: List[str]):
        """
        Serialize DefectComplex to graph.

        Parameters
        ----------
        graph : KnowledgeGraph
            The knowledge graph instance
        sample_id : str
            The sample URI
        defect_uris : List[str]
            List of actual defect URIs collected during serialization
        """
        name = sample_id
        material = get_material(graph, sample_id)
        defect_complex = graph.create_node(
            f"{name}_DefectComplex", UNSAFECDCO.DefectComplex
        )
        graph.add((material, UNSAFECDCO.hasDefectComplex, defect_complex))

        # Link the actual defect URIs
        if defect_uris:
            for defect_uri in defect_uris:
                graph.add(
                    (
                        URIRef(defect_uri),
                        UNSAFECDCO.isPartOfDefectComplex,
                        defect_complex,
                    )
                )

        # Add relative_distance if present
        if self.relative_distance is not None:
            graph.add(
                (
                    defect_complex,
                    UNSAFECDCO.hasRelativeDistance,
                    Literal(self.relative_distance),
                ),
                validate=False,
            )

    @classmethod
    def from_graph(cls, graph, sample):
        """
        Deserialize DefectComplex from graph.

        Returns
        -------
        DefectComplex or None
            DefectComplex instance if found in graph, None otherwise
        """
        material = get_material(graph, sample)

        # Find the DefectComplex node
        defect_complex_node = graph.value(material, UNSAFECDCO.hasDefectComplex)

        if defect_complex_node is None:
            return None

        # Get relative distance if present
        relative_distance = graph.value(
            defect_complex_node, UNSAFECDCO.hasRelativeDistance
        )
        relative_distance_str = (
            relative_distance.toPython() if relative_distance else None
        )

        # Find all defects that are part of this complex
        defect_types = []
        for defect_uri in graph.subjects(
            UNSAFECDCO.isPartOfDefectComplex, defect_complex_node
        ):
            # Extract defect type from URI (e.g., "sample:abc_StackingFault" -> "stacking_fault")
            defect_uri_str = str(defect_uri)

            # Map class names back to field names
            class_to_field = {
                "PointDefect": "point_defect",
                "Vacancy": "vacancy",
                "Substitutional": "substitutional",
                "Interstitial": "interstitial",
                "Dislocation": "dislocation",
                "EdgeDislocation": "edge_dislocation",
                "ScrewDislocation": "screw_dislocation",
                "MixedDislocation": "mixed_dislocation",
                "StackingFault": "stacking_fault",
                "GrainBoundary": "grain_boundary",
                "TiltGrainBoundary": "tilt_grain_boundary",
                "TwistGrainBoundary": "twist_grain_boundary",
                "SymmetricTiltGrainBoundary": "symmetric_tilt_grain_boundary",
                "MixedGrainBoundary": "mixed_grain_boundary",
            }

            # Extract class name from URI
            for class_name, field_name in class_to_field.items():
                if defect_uri_str.endswith(f"_{class_name}"):
                    defect_types.append(field_name)
                    break

        return cls(ids=defect_types, relative_distance=relative_distance_str)
