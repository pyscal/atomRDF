"""
https://docs.python.org/3/library/operator.html
"""
from rdflib import Namespace
import numbers

def strip_name(uri, delimiter, get_what='name'):
    if get_what == "name": 
        if delimiter == "/":
            uri_split = uri.split(delimiter)
            if len(uri_split)>1:
                name = ":".join(uri_split[-2:])
            else:
                name = uri 
        else:
            uri_split = uri.split(delimiter)
            name = uri_split[-1]
            uri_split = uri_split[0].split("/")
            if len(uri_split)>0:
                namespace = uri_split[-1]
                name = ":".join([namespace, name])
        return name
    
    elif get_what == "namespace":
        if delimiter == "/":
            uri_split = uri.split(delimiter)
            if len(uri_split)>1:
                namespace = uri_split[-2]
            else:
                namespace = uri 
        else:
            uri_split = uri.split(delimiter)
            uri_split = uri_split[0].split("/")
            if len(uri_split)>0:
                namespace = uri_split[-1]
            else:
                namespace = uri
        return namespace 


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
        self._condition = None

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
        return strip_name(self.uri, self.delimiter, get_what="name")        

    @property
    def namespace(self):
        return strip_name(self.uri, self.delimiter, get_what="namespace")        

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

    @property
    def query_name(self):
        """
        What it is called in a sparql query
        """
        if self.node_type == "data_property":
            return self.name + "value"
        return self.name

    @property
    def query_name_without_prefix(self):
        """
        What it is called in a sparql query
        """
        if self.node_type == "data_property":
            return self.name_without_prefix + "value"
        return self.name_without_prefix

    def __repr__(self):
        return str(self.name)

    #convenience methods for overload checking
    def _ensure_condition_exists(self):
        if self._condition is None:
            raise ValueError("Individual terms should have condition for this operation!")

    def _is_term(self, val):
        if not isinstance(val, OntoTerm):
            raise TypeError("can only be performed with an OntoTerm!")

    def _is_number(self, val):
        if not isinstance(val, numbers.Number):
            raise TypeError("can only be performed with a number!")
    
    def _is_data_node(self):
        if not self.node_type == "data_property":
            raise TypeError("This operation can only be performed with a data property!")

    def _create_condition_string(self, condition, val):
        return f'(?{self.query_name_without_prefix}{condition}\"{val}\"^^xsd:{self.range[0]})'
    
    #overloading operators
    def __eq__(self, val):
        """
        = 
        """
        #self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string("=", val)
        return self

    def __lt__(self, val):
        self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string("<", val)
        return self

    def __le__(self, val):
        self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string("<=", val)
        return self


    def __ne__(self, val):
        #self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string("!=", val)
        return self
    
    def __ge__(self, val):
        self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string(">=", val)
        return self
    
    def __gt__(self, val):
        self._is_number(val)
        self._is_data_node()
        self._condition = self._create_condition_string(">", val)
        return self

    def __and__(self, term):
        self._is_term(term)
        self._is_data_node()
        term._is_data_node()
        self._ensure_condition_exists()
        term._ensure_condition_exists()        
        self._condition = "&&".join([self._condition, term._condition])
        self._condition = f'({self._condition})'
        return self

    def and_(self, term):
        self.__and__(term)

    def __or__(self, term):
        self._is_term(term)
        self._is_data_node()
        term._is_data_node()
        self._ensure_condition_exists()
        term._ensure_condition_exists()        
        self._condition = "||".join([self._condition, term._condition])
        self._condition = f'({self._condition})'
        return self

    def or_(self, term):
        self.__or__(term)