"""
https://docs.python.org/3/library/operator.html
"""

from rdflib import Namespace, URIRef
import numbers
import copy

def _get_namespace_and_name(uri):
    uri_split = uri.split('#')
    if len(uri_split) > 1:
        #possible that delimiter is #
        name = uri_split[-1]
        namespace_split = uri_split[0].split('/')
        namespace = namespace_split[-1]
    else:
        uri_split = uri.split('/')
        if len(uri_split) > 1:
            name = uri_split[-1]
            namespace = uri_split[-2]
        else:
            name = uri_split[-1]
            namespace = ""
    return namespace, name

def _get_namespace_with_prefix(uri):
    uri_split = uri.split('#')
    if len(uri_split) > 1:
        #possible that delimiter is #
        namespace = uri_split[0]
    else:
        uri_split = uri.split('/')
        if len(uri_split) > 1:
            namespace = "/".join(uri_split[:-1])
        else:
            namespace = ""
        if namespace[-1] != "#":
            namespace += "/"
    return namespace

def strip_name(uri, get_what="name", namespace=None):
    if namespace is None:
        namespace, name = _get_namespace_and_name(uri)
    else:
        _, name = _get_namespace_and_name(uri)
    if get_what == "namespace":
        return namespace

    elif get_what == "name":
        return ":".join([namespace, name])


class OntoTerm:
    def __init__(
        self,
        uri = None,
        namespace=None,
        node_type=None,
        dm=[],
        rn=[],
        data_type=None,
        node_id=None,
        delimiter="/",
        description=None,
        label=None,
    ):

        self.uri = uri
        # type: can be object property, data property, or class
        self.node_type = node_type
        # now we need domain and range
        self.domain = dm
        self.range = rn
        # datatype, that is only need for data properties
        self.data_type = data_type
        # identifier
        self.node_id = node_id
        self.associated_data_node = None
        self.subclasses = []
        self.named_individuals = []
        self.equivalent_classes = []
        self.subproperties = []
        self.delimiter = delimiter
        self.description = description
        self.label = label
        self._description = None
        self._label = None
        self.is_domain_of = []
        self.is_range_of = []
        self._condition = None
        self._namespace = namespace
        # name of the class
        self._name = None
        #parents for the class; these are accumulated
        #when using the >> operator
        self._parents = []
        #condition parents are the parents that have conditions
        #these are accumulated when using the & or || operators
        self._condition_parents = []
        self._object = None
    
    @property
    def URIRef(self):
        return URIRef(self._uri)

    @property
    def uri(self):
        """
        Get the URI of the ontology term.

        Returns
        -------
        str
            The URI of the ontology term.
        """
        return self._uri

    @uri.setter
    def uri(self, val):
        self._uri = val

    @property
    def description(self):
        """
        Get the description of the term.

        Returns
        -------
        str
            The description of the term.
        """
        return self._description
    
    @description.setter
    def description(self, val):
        if isinstance(val, list):
            if len(val) == 0:
                val = ''
            elif len(val) > 1:
                val = ". ".join(val)
            else:
                val = val[0]
        self._description = val

    @property
    def label(self):
        """
        Get the label of the term.

        Returns
        -------
        str
            The label of the term.
        """
        return self._label

    @label.setter
    def label(self, val):
        if isinstance(val, list):
            if len(val) == 0:
                val = ''
            elif len(val) > 1:
                val = ". ".join(val)
            else:
                val = val[0]
        self._label = val

    @property
    def name(self):
        """
        Get the name of the term.

        Returns
        -------
        str
            The name of the term.
        """
        if self._name is not None:
            return self._name
        return strip_name(
            self.uri, get_what="name",namespace=self.namespace,
        )
    
    @name.setter
    def name(self, val):
        self._name = val

    @property
    def name_without_prefix(self):
        """
        Get the name without the namespace prefix.

        Returns
        -------
        str
            The name of the term without the namespace prefix.
        """
        name = self.name
        name = name.replace("–", "")
        name = name.replace("-", "")
        name = name.split(':')[-1]
        return name

    @property
    def namespace(self):
        """
        Get the namespace of the term.

        Returns
        -------
        str
            The namespace of the term.
        """
        if self._namespace is not None:
            return self._namespace
        else:
            return strip_name(self.uri, get_what="namespace")

    @property
    def namespace_with_prefix(self):
        """
        Get the namespace of the term with the prefix.
    
        Returns
        -------
        str
            The namespace of the term with the prefix.
        """
        return _get_namespace_with_prefix(self.uri)

    @property
    def namespace_object(self):
        """
        Get the namespace object for the term.

        Returns
        -------
        object
            The namespace object for the term.

        """
        return self.URIRef

    @property
    def query_name(self):
        """
        Get the name of the term as it appears in a SPARQL query.

        Returns
        -------
        str
            The name of the term in a SPARQL query.

        Notes
        -----
        If the term is a data property, the name will be appended with "value".

        """
        if self.node_type == "data_property":
            return self.name + "value"
        elif self.node_type == "object_property":
            if len(self.range) > 0:
                #this has a domain
                return self.range[0]
        return self.name

    @property
    def variable_name(self):
        """
        Get the name of the term to use as a variable in a SPARQL query.

        Returns
        -------
        str
            The name of the term in a SPARQL query.

        """
        name_list = [x.name_without_prefix for x in self._parents]
        name_list.append(self.name_without_prefix)
        name =  "_".join(name_list) 
        
        if self.node_type == "data_property":
            return name + "value"
        
        return name

    @property
    def query_name_without_prefix(self):
        """
        Get the name of the term as it appears in a SPARQL query without prefix.

        Returns
        -------
        str
            The name of the term in a SPARQL query.

        Notes
        -----
        If the term is a data property, the name will be suffixed with "value".
        """
        if self.node_type == "data_property":
            return self.name_without_prefix + "value"
        return self.name_without_prefix

    def toPython(self):
        return self.uri

    def __repr__(self):
        if self.description is not None:
            return str(self.name + '\n' + self.description)
        else:
            return str(self.name)

    def _clean_datatype(self, r):
        if r == "str":
            return "string"
        return r

    # convenience methods for overload checking
    def _ensure_condition_exists(self):
        if self._condition is None:
            raise ValueError(
                "Individual terms should have condition for this operation!"
            )

    def _is_term(self, val):
        if not isinstance(val, OntoTerm):
            raise TypeError("can only be performed with an OntoTerm!")

    def _is_number(self, val):
        if not isinstance(val, numbers.Number):
            raise TypeError("can only be performed with a number!")

    def _is_data_node(self):
        if not self.node_type == "data_property":
            raise TypeError(
                "This operation can only be performed with a data property!"
            )
    


    def _create_condition_string(self, condition, val):
        return f'(?{self.variable_name}{condition}"{val}"^^xsd:{self._clean_datatype(self.range[0])})'

    # overloading operators
    def __eq__(self, val):
        """
        =
        """
        #print("eq")
        #print(f'lhs {self} rhs {val}')
        self._is_data_node()
        item = copy.deepcopy(self)
        item._condition = item._create_condition_string("=", val)
        return item

    def __lt__(self, val):
        self._is_number(val)
        self._is_data_node()
        item = copy.deepcopy(self)
        item._condition = item._create_condition_string("<", val)
        return item

    def __le__(self, val):
        self._is_number(val)
        self._is_data_node()
        item = copy.deepcopy(self)
        item._condition = item._create_condition_string("<=", val)
        return item

    def __ne__(self, val):
        # self._is_number(val)
        self._is_data_node()
        item = copy.deepcopy(self)
        item._condition = item._create_condition_string("!=", val)
        return item

    def __ge__(self, val):
        self._is_number(val)
        self._is_data_node()
        item = copy.deepcopy(self)
        item._condition = item._create_condition_string(">=", val)
        return item

    def __gt__(self, val):
        #print("gt")
        #print(f'lhs {self} rhs {val}')
        self._is_number(val)
        self._is_data_node()
        item = copy.deepcopy(self)
        self._condition = self._create_condition_string(">", val)
        return self

    def __and__(self, term):
        self._is_term(term)
        self._is_data_node()
        term._is_data_node()
        self._ensure_condition_exists()
        term._ensure_condition_exists()
        item = copy.deepcopy(self)
        item._condition = "&&".join([item._condition, term._condition])
        item._condition = f"({item._condition})"
        item._condition_parents.append(copy.deepcopy(term))
        #and clean up the inbound term
        if item.name != term.name:
            term.refresh_condition()
        return item

    def and_(self, term):
        self.__and__(term)

    def __or__(self, term):
        self._is_term(term)
        self._is_data_node()
        term._is_data_node()
        self._ensure_condition_exists()
        term._ensure_condition_exists()
        item = copy.deepcopy(self)
        item._condition = "||".join([item._condition, term._condition])
        item._condition = f"({item._condition})"
        item._condition_parents.append(copy.deepcopy(term))
        #and clean up the inbound term
        if item.name != term.name:
            term.refresh_condition()
        return item

    def or_(self, term):
        self.__or__(term)

    def __matmul__(self, term):
        #print("matmul")
        #print(f'lhs {self} rhs {term}')
        item = copy.deepcopy(self)
        item._parents.append(copy.deepcopy(term))
        return item
    
    def refresh_condition(self):
        self._condition = None
        self._condition_parents = []

    def refresh(self):
        self._condition = None
        self._parents = []
        self._condition_parents = []

