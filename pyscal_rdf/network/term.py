from rdflib import Namespace

class OntoTerm:
    def __init__(self, uri, node_type=None, 
                dm=[], rn=[], data_type=None, 
                 node_id=None, delimiter='/'):
        """
        This is class that represents an ontology element
        """
        self.uri = uri
        #name of the class
        self._name = None
        #type: can be object property, data property, or class
        self.node_type = node_type
        #now we need domain and range
        self.domain = dm
        self.range = rn
        #datatype, that is only need for data properties
        self.data_type = data_type
        #identifier
        self.node_id = node_id
        self.subclasses = []
        self.delimiter = delimiter
        self.is_domain_of = []
        self.is_range_of = []

    @property
    def uri(self):
        return self._uri
    
    @uri.setter
    def uri(self, val):
        self._uri = val
    
    @property
    def name_without_prefix(self):
        uri_split = self.uri.split(self.delimiter)
        if len(uri_split)>0:
            return uri_split[-1]
        else:
            return self.uri

    @property
    def name(self):
        uri_split = self.uri.split(self.delimiter)
        if len(uri_split)>1:
            return ":".join(uri_split[-2:])
        else:
            return self.uri
    
    @property
    def namespace(self):
        uri_split = self.uri.split(self.delimiter)
        if len(uri_split)>1:
            return uri_split[-2]
        else:
            return self.uri

    @property
    def namespace_with_prefix(self):
        uri_split = self.uri.split(self.delimiter)
        if len(uri_split)>1:
            namespace = self.delimiter.join(uri_split[:-1]) + self.delimiter
            return namespace
        else:
            return self.uri
                    
    @property
    def namespace_object(self):
        uri_split = self.uri.split(self.delimiter)
        if len(uri_split)>1:
            namespace = self.delimiter.join(uri_split[:-1]) + self.delimiter
            prop = uri_split[-1]
            return getattr(Namespace(namespace), prop)
        else:
            return self.uri

    def __repr__(self):
        return str(self.name)