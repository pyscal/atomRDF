import networkx as nx
import graphviz
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
from pyscal_rdf.network.parser import OntoParser
from pyscal_rdf.network.term import OntoTerm, strip_name
from pyscal3.atoms import AttrSetter

owlfile = os.path.join(os.path.dirname(__file__), "../data/cmso.owl")

def _replace_name(name):
    return ".".join(name.split(':'))

class OntologyNetwork:
    """
    Network representation of Onto
    """
    def __init__(self, infile=None, delimiter='/'):
        if infile is None:
            infile = owlfile
            
        self.g = nx.DiGraph()
        self.onto = OntoParser(infile, delimiter=delimiter)
        self.onto.attributes['data_node'] = []
        self.data_prefix = 'value'
        self.terms = AttrSetter()
        self._parse_all()        
    
    def _assign_attributes(self):
        mapdict = {}
        #add first level - namespaces
        for key in self.namespaces.keys():
            mapdict[key] = {}
        
        #now iterate over all attributes
        for k1 in ['class', 'object_property', 'data_property']:
            for k2, val in self.onto.attributes[k1].items():
                mapdict[val.namespace][val.name_without_prefix] = val
        

        self.terms._add_attribute(mapdict) 

    def _parse_all(self):
        #call methods
        self._add_class_nodes()
        self._add_object_properties()
        self._add_data_properties()
        self._assign_attributes()

    def __add__(self, ontonetwork):
        #add onto network
        self.onto = self.onto + ontonetwork.onto
        #now parse again
        self._parse_all()
        return self

    def strip_name(self, name):
        raw = name.split(':')
        if len(raw) > 1:
            return raw[-1]
        return name

    @property
    def attributes(self):
        return self.onto.attributes

    @property
    def namespaces(self):
        return self.onto.namespaces

    @property
    def extra_namespaces(self):
        return self.onto.extra_namespaces

    def __radd__(self, ontonetwork):
        return self.__add__(ontonetwork)

    def get_shortest_path(self, source, target, triples=False):
        path = nx.shortest_path(self.g, source=source, target=target)
        if triples:
            triple_list = []
            for x in range(len(path)//2):
                triple_list.append(path[2*x:2*x+3])
            return triple_list
        return path
    
    def _add_class_nodes(self):
        for key, val in self.onto.attributes['class'].items():
            self.g.add_node(val.name, node_type='class')
    
    def _add_object_properties(self):
        for key, val in self.onto.attributes['object_property'].items():
            self.g.add_node(val.name, node_type='object_property')
            #find domain
            for d in val.domain:
                self.g.add_edge(d, val.name)
            for r in val.range:
                self.g.add_edge(val.name, r)

    
    def _add_data_properties(self):
        for key, val in self.onto.attributes['data_property'].items():
            self.g.add_node(val.name, node_type='data_property')
            for d in val.domain:
                self.g.add_edge(d, val.name)
            for r in val.range:
                data_node = f'{val.name}{self.data_prefix}'
                self.onto.attributes['data_node'].append(data_node)
                self.g.add_node(data_node, node_type='literal', data_type=r)
                self.g.add_edge(val.name, data_node)

                
    def add_namespace(self, namespace_name, namespace_iri):
        """
        Add a new namespace
        """
        if namespace_name not in self.onto.namespaces.keys():
            self.onto.namespaces[namespace_name] = namespace_iri
        else:
            raise KeyError("namespace is already there!")

    def add_term(self, uri, node_type, namespace=None,
                dm=[], rn=[], data_type=None, 
                node_id=None, delimiter='/'):
        """
        Add a node
        """
        #namespace = strip_name(uri, delimiter, get_what="namespace")
        #name = strip_name(uri, delimiter, get_what="name")
        term = OntoTerm(uri, namespace=namespace,
            node_type=node_type, dm =dm, 
            rn=rn, data_type=data_type, node_id=node_id,
            delimiter=delimiter)
        if not term.namespace in self.onto.namespaces.keys():
            raise ValueError("Namespace not found, first add namespace")
        self.onto.attributes[node_type][term.name] = term
        self._assign_attributes()

    def add_path(self, triple):
        """
        Add a triple as path. Note that all attributes of the triple should already
        exist in the graph. The ontology itself is not modified. Only the graph
        representation of it is.
        The expected use is to bridge between two(or more) different ontologies.
        Therefore, mapping can only be between classes.
        """
        sub = triple[0]
        pred = triple[1]
        obj = triple[2]

        if sub not in self.onto.attributes['class'].keys():
            raise ValueError(f'{sub} not found in self.attributes')

        #now add
        subclasses = self.onto._get_subclasses(sub)
        for subclass in subclasses:
            self.g.add_edge(subclass, pred)            
        
        #now add pred
        if pred in self.onto.attributes['object_property'].keys():
            if obj not in self.onto.attributes['class'].keys():
                raise ValueError(f'{obj} not found in self.attributes')
            subclasses = self.onto._get_subclasses(obj)
            for subclass in subclasses:
                self.g.add_edge(pred, subclass) 

        #another possibility it is data property
        elif pred in self.onto.attributes['data_property'].keys():
            data_node = f'{pred}{self.data_prefix}'
            self.g.add_node(data_node, node_type='literal', data_type=obj)
            self.g.add_edge(pred, data_node)
        
        else:
            raise ValueError(f'{pred} not found in self.attributes')

    def draw(self, styledict = {"class": {"shape":"box"},
                                "object_property": {"shape":"ellipse"},
                                "data_property": {"shape":"ellipse"},
                                "literal": {"shape":"parallelogram"},}):
        dot = graphviz.Digraph()
        node_list = list(self.g.nodes(data='node_type'))
        edge_list = list(self.g.edges)
        for node in node_list:
            name = _replace_name(node[0])
            if node[1] is not None:
                t = node[1]
                dot.node(name, shape=styledict[t]['shape'], fontsize="6")
        for edge in edge_list:
            dot.edge(_replace_name(edge[0]), _replace_name(edge[1]))
        return dot

    def get_path_from_sample(self, target):
        path = self.get_shortest_path(source="cmso:ComputationalSample", target=target, triples=True)
        return path
        
        
    def create_query(self, source, destinations, condition=None, enforce_types=True):
        """
        values is a dict with keys value, operation
        """
        if not isinstance(destinations, list):
            destinations = [destinations]
        
        source_name = source.query_name
        destination_names = [destination.query_name for destination in destinations]
        
        #if condition is specified, and is not there, add it
        if condition is not None:
            if condition.query_name not in destination_names:
                destination_names.append(condition.query_name)
        
        #add source if not available
        if source_name not in destination_names:
            destination_names.append(source_name)

        #start prefix of query
        query = []
        for key, val in self.namespaces.items():
            query.append(f'PREFIX {key}: <{val}>')
        for key, val in self.extra_namespaces.items():
            query.append(f'PREFIX {key}: <{val}>')
        
        #now for each destination, start adding the paths in the query
        all_triplets = {}
        for destination in destination_names:
            triplets = self.get_shortest_path(source_name, destination, triples=True)
            all_triplets[destination] = triplets
        
        select_destinations = [f'?{self.strip_name(destination)}' for destination in destination_names]
        query.append(f'SELECT DISTINCT {" ".join(select_destinations)}')
        query.append("WHERE {")
        
        #now add corresponding triples
        for destination in destination_names:
            for triple in all_triplets[destination]:
                #print(triple)
                query.append("    ?%s %s ?%s ."%(self.strip_name(triple[0]), 
                                                 triple[1], 
                                                 self.strip_name(triple[2])))
        
        #we enforce types of the source and destination
        if enforce_types:
            if source.node_type == 'class':
                query.append("    ?%s rdf:type %s ."%(self.strip_name(source.query_name), source.query_name))
            for destination in destinations:
                if destination.node_type == 'class':
                    query.append("    ?%s rdf:type %s ."%(self.strip_name(destination.query_name), destination.query_name))
        #now we have to add filters
        #filters are only needed if it is a dataproperty
        filter_text = ''
        
        #make filters; get all the unique filters from all the classes in destinations
        if condition is not None:
            if condition._condition is not None:
                filter_text = condition._condition

        if filter_text != '':
            query.append(f'FILTER {filter_text}')
        query.append('}')
        return '\n'.join(query)

        
            
            
            

        
        