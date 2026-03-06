"""
This module provides the Namespace class for managing namespaces in the AtomRDF library.

The Namespace class extends the rdflib.Namespace class and provides additional functionality for working with namespaces.

Classes
-------
Namespace
    A class representing a namespace in the AtomRDF library.
"""

import os
import json
import numpy as np
from rdflib import URIRef
from rdflib import Namespace as RDFLibNamespace
from rdflib import Literal as RDFLibLiteral
from pyscal3.atoms import AttrSetter
from tools4rdf.network.network import OntologyNetwork


def Literal(value, datatype=None):
    if datatype is not None:
        return RDFLibLiteral(value, datatype=datatype)
    elif isinstance(value, list):
        return RDFLibLiteral(json.dumps(value))
    elif isinstance(value, np.ndarray):
        return RDFLibLiteral(json.dumps(value.tolist()))
    else:
        return RDFLibLiteral(value)


class Namespace(AttrSetter, RDFLibNamespace):
    """A class representing a namespace in the AttrSetter library.

    This class extends the `rdflib.Namespace` classes.

    Parameters
    ----------
    infile : str
        The input file path for the OWL ontology.
    delimiter : str, optional
        The delimiter used in the input file. Defaults to "/".

    Attributes
    ----------
    network : OntologyNetwork
        The ontology network associated with the namespace.
    name : str
        The name of the namespace.
    """

    def __new__(cls, infile, delimiter="/"):
        """
        Create the Namespace with the canonical base IRI from the OWL file.

        ``rdflib.Namespace`` (and ``rdflib.URIRef``) inherit from ``str``, so
        the actual IRI string is immutable after object creation — it must be
        set here in ``__new__``, not in ``__init__``.  Without this, the
        namespace IRI was the local file *path* instead of the purls HTTP URI.
        """
        network = OntologyNetwork(infile)
        base_iri = network.onto.base_iri
        # Normalise: rdflib.Namespace builds term URIs as base_iri + term_name,
        # so a trailing '/' (or '#') is required for correct term construction.
        if not base_iri.endswith("/") and not base_iri.endswith("#"):
            base_iri = base_iri + "/"
        instance = super().__new__(cls, base_iri)
        # Cache the loaded network to avoid parsing the OWL file twice.
        instance._loaded_network = network
        return instance

    def __init__(self, infile, delimiter="/"):
        """
        Initialize the Namespace class.

        Parameters
        ----------
        infile : str
            The input file path for the OWL ontology.
        delimiter : str, optional
            The delimiter used in the input file. Defaults to "/".
        """
        AttrSetter.__init__(self)
        # Reuse the network loaded in __new__ to avoid parsing twice.
        self.network = self._loaded_network
        self.name = self.network.onto.base_iri.split("/")[-1]
        mapdict = {}

        # now iterate over all attributes
        for k1 in ["class", "object_property", "data_property"]:
            for k2, val in self.network.onto.attributes[k1].items():
                mapdict[val.name_without_prefix] = val

        # add attributes
        self._add_attribute(mapdict)


file_location = os.path.dirname(__file__)

CMSO = Namespace(os.path.join(file_location, "data/cmso.owl"))
LDO = Namespace(os.path.join(file_location, "data/ldo.owl"))
PLDO = Namespace(os.path.join(file_location, "data/pldo.owl"))
PODO = Namespace(os.path.join(file_location, "data/podo.owl"))
ASMO = Namespace(os.path.join(file_location, "data/asmo.owl"))
MATH = Namespace(os.path.join(file_location, "data/asmo.owl"))
CDCO = Namespace(os.path.join(file_location, "data/cdco.owl"))

PROV = RDFLibNamespace("http://www.w3.org/ns/prov#")
MDO = RDFLibNamespace("https://w3id.org/mdo/calculation/")
DCAT = RDFLibNamespace("http://www.w3.org/ns/dcat#")

# Plain RDFLibNamespace aliases for when you need simple string-concatenation
# access without OWL attribute validation (e.g. term names that differ from the
# OWL class name).  Now use the canonical purls IRIs from the fixed Namespace
# objects above instead of hard-coded strings.
UNSAFECMSO = RDFLibNamespace(str(CMSO))
UNSAFEASMO = RDFLibNamespace(str(ASMO))
UNSAFECDCO = RDFLibNamespace(str(CDCO))
