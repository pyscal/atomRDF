import os
from rdflib import Literal, URIRef
from rdflib import Namespace as RDFLibNamespace
from pyscal3.atoms import AttrSetter

from atomrdf.network.network import OntologyNetwork

class Namespace(AttrSetter, RDFLibNamespace):
    def __init__(self, infile, delimiter='/'):
        AttrSetter.__init__(self)
        self.network = OntologyNetwork(infile, delimiter=delimiter)
        #print(type(self.network.onto.tree.base_iri))
        #self.namespace = RDFLibNamespace(self.network.onto.tree.base_iri)
        RDFLibNamespace.__init__(self.network.onto.tree.base_iri)
        #self.namespace = RDFLibNamespace("http://purls.helmholtz-metadaten.de/cmso/")
        self.name = self.network.onto.tree.name
        mapdict = {}
        
        #now iterate over all attributes
        for k1 in ['class', 'object_property', 'data_property']:
            for k2, val in self.network.onto.attributes[k1].items():
                if val.namespace == self.name:
                    mapdict[val.name_without_prefix] = val
        
        #add attributes
        self._add_attribute(mapdict)


file_location = os.path.dirname(__file__)

CMSO = Namespace(os.path.join(file_location,  'data/cmso.owl'))
PLDO = Namespace(os.path.join(file_location,  'data/pldo.owl'))
PODO = Namespace(os.path.join(file_location,  'data/podo.owl'))
ASMO = Namespace(os.path.join(file_location,  'data/asmo.owl'))
PROV = RDFLibNamespace("http://www.w3.org/ns/prov#")