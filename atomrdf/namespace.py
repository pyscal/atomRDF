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

from atomrdf.network.network import OntologyNetwork

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
    """A class representing a namespace in the AtomRDF library.

    This class extends the `rdflib.Namespace` classes.

    Parameters
    ----------
    infile : str
        The input file path.
    delimiter : str, optional
        The delimiter used in the input file. Defaults to "/".

    Attributes
    ----------
    network : OntologyNetwork
        The ontology network associated with the namespace.
    name : str
        The name of the namespace.
    """

    def __init__(self, infile, delimiter="/"):
        """
        Initialize the Namespace class.

        Parameters
        ----------
        infile : str
            The input file path.
        delimiter : str, optional
            The delimiter used in the input file. Defaults to "/".
        """
        AttrSetter.__init__(self)
        self.network = OntologyNetwork(infile, delimiter=delimiter)
        RDFLibNamespace.__init__(self.network.onto.tree.base_iri)
        self.name = self.network.onto.tree.name
        mapdict = {}

        # now iterate over all attributes
        for k1 in ["class", "object_property", "data_property"]:
            for k2, val in self.network.onto.attributes[k1].items():
                if val.namespace == self.name:
                    mapdict[val.name_without_prefix] = val

        # add attributes
        self._add_attribute(mapdict)


file_location = os.path.dirname(__file__)

CMSO = Namespace(os.path.join(file_location, "data/cmso.owl"))
LDO = Namespace(os.path.join(file_location, "data/ldo.owl"))
PLDO = Namespace(os.path.join(file_location, "data/pldo.owl"))
PODO = Namespace(os.path.join(file_location, "data/podo.owl"))
ASMO = Namespace(os.path.join(file_location, "data/asmo.owl"))
PROV = RDFLibNamespace("http://www.w3.org/ns/prov#")
MDO = RDFLibNamespace("https://w3id.org/mdo/calculation/")

MATH = RDFLibNamespace("http://purls.helmholtz-metadaten.de/asmo/")
UNSAFECMSO = RDFLibNamespace("http://purls.helmholtz-metadaten.de/cmso/")
UNSAFEASMO = RDFLibNamespace("http://purls.helmholtz-metadaten.de/asmo/")
