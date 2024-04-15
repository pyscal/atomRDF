from rdflib import Namespace, Literal, URIRef
from pyscal3.atoms import AttrSetter

from atomrdf.network.network import OntologyNetwork

class StrictNamespace(AttrSetter):
    def __init__(self, infile, delimiter='/'):
        AttrSetter.__init__(self)
        self.network = OntologyNetwork(infile, delimiter=delimiter)
        self.namespace = Namespace(self.network.onto.tree.base_iri)
        self.name = self.network.onto.tree.name
        mapdict = {}
        
        #now iterate over all attributes
        for k1 in ['class', 'object_property', 'data_property']:
            for k2, val in self.network.onto.attributes[k1].items():
                if val.namespace == self.name:
                    mapdict[val.name_without_prefix] = val
        
        #add attributes
        self._add_attribute(mapdict)