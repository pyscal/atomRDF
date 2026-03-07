from typing import List, Optional
import uuid
from pydantic import Field
from rdflib import URIRef

from atomrdf.datamodels.basemodels import BaseModel, TemplateMixin
from atomrdf.namespace import DCAT, DCTERMS, FOAF, Literal


class Creator(BaseModel, TemplateMixin):
    """A person who created a dataset, mapped to foaf:Person.

    ``id`` (from TemplateMixin): URI identifying the person, e.g. an ORCID.
    """

    name: Optional[str] = Field(default=None, description="Full name (foaf:name)")

    def to_graph(self, graph) -> URIRef:
        uri = self.id if self.id else f"person:{uuid.uuid4()}"
        subject = graph.create_node(uri, FOAF.Person, label=self.name)
        if self.name:
            graph.add((subject, FOAF.name, Literal(self.name)))
        return subject

    @classmethod
    def from_graph(cls, graph, subject: URIRef) -> "Creator":
        name_node = graph.value(subject, FOAF.name)
        return cls(
            id=str(subject),
            name=str(name_node) if name_node is not None else None,
        )


class Publication(BaseModel, TemplateMixin):
    """A bibliographic resource (paper), mapped to dcterms:BibliographicResource.

    ``id`` (from TemplateMixin): URI identifying the paper.
    """

    identifier: Optional[str] = Field(default=None, description="DOI string")
    title: Optional[str] = Field(default=None, description="dcterms:title")

    def to_graph(self, graph) -> URIRef:
        uri = self.id if self.id else f"publication:{uuid.uuid4()}"
        subject = graph.create_node(
            uri, DCTERMS.BibliographicResource, label=self.title
        )
        if self.identifier:
            graph.add((subject, DCTERMS.identifier, Literal(self.identifier)))
        if self.title:
            graph.add((subject, DCTERMS.title, Literal(self.title)))
        return subject

    @classmethod
    def from_graph(cls, graph, subject: URIRef) -> "Publication":
        identifier_node = graph.value(subject, DCTERMS.identifier)
        title_node = graph.value(subject, DCTERMS.title)
        return cls(
            id=str(subject),
            identifier=str(identifier_node) if identifier_node is not None else None,
            title=str(title_node) if title_node is not None else None,
        )


class Dataset(BaseModel, TemplateMixin):
    """A dataset, mapped to dcat:Dataset."""

    identifier: Optional[str] = Field(
        default=None, description="DOI URI for the dataset"
    )
    title: Optional[str] = Field(default=None, description="dcterms:title")
    creators: List[Creator] = Field(default_factory=list)
    publication: Optional[Publication] = Field(default=None)
    samples: List[str] = Field(
        default_factory=list, description="Sample URIs linked via dcterms:isPartOf"
    )

    def to_graph(self, graph) -> URIRef:
        uri = self.identifier if self.identifier else f"dataset:{uuid.uuid4()}"
        subject = graph.create_node(uri, DCAT.Dataset, label=self.title)
        if self.identifier:
            graph.add((subject, DCTERMS.identifier, Literal(self.identifier)))
        if self.title:
            graph.add((subject, DCTERMS.title, Literal(self.title)))
        for creator in self.creators:
            person_uri = creator.to_graph(graph)
            graph.add((subject, DCTERMS.creator, person_uri))
        if self.publication is not None:
            pub_uri = self.publication.to_graph(graph)
            graph.add((subject, DCTERMS.isReferencedBy, pub_uri))
        for sample_uri in self.samples:
            graph.add((URIRef(sample_uri), DCTERMS.isPartOf, subject))
        return subject

    @classmethod
    def from_dict(cls, data: dict) -> "Dataset":
        """Construct a Dataset from a plain dictionary (e.g., parsed YAML)."""
        creators = [Creator(**c) for c in data.get("creators", [])]
        pub_data = data.get("publication")
        publication = Publication(**pub_data) if pub_data else None
        return cls(
            identifier=data.get("identifier"),
            title=data.get("title"),
            creators=creators,
            publication=publication,
            samples=data.get("samples", []),
        )
