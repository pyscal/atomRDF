import networkx as nx
import graphviz
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
from pyscal_rdf.network.parser import OntoParser

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
        self._parse_all()
        
        
    def _parse_all(self):
        #call methods
        self._add_class_nodes()
        self._add_object_properties()
        self._add_data_properties()

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
        
    def phrase_to_sparql(self, phrase):
        def _extract_operation(phr):
            r = phr.split(' ')
            if len(r) != 3:
                raise ValueError('wrong filters!')
            return f'?value{r[1]}\"{r[2]}\"^^xsd:datatype'

        conditions = []
        operation = None
        
        raw = phrase.split(' and ')
        
        if len(raw) > 1:
            operation = '&&'
        if operation is None:
            raw = phrase.split(' or ')
            if len(raw) > 1:
                operation = '||'
        
        if operation is not None:
            for ph in raw:
                conditions.append(_extract_operation(ph))
        else:
            conditions.append(_extract_operation(phrase))
        full_str = f' {operation} '.join(conditions)
        #replace values
        return full_str
        

    def validate_values(self, destinations, values):
        combinator_dict = {'and': '&&', 'or': '||'}
        combinator_list = values[1::2]
        phrase_list = values[::2]
        if not len(combinator_list) == len(destinations)-1:
            raise ValueError("Invalid combinations!")
        
        sparql_phrase_list = []
        for phrase, destination in zip(phrase_list, destinations):
            sparql_phrase = self.phrase_to_sparql(phrase)
            sparql_phrase = sparql_phrase.replace('value', self.strip_name(destination))
            sparql_phrase = sparql_phrase.replace('datatype', self.g.nodes[destination]['data_type'])
            sparql_phrase_list.append(sparql_phrase)
            
        #combine phrases with phrase list
        updated_sparql_phrase_list = []
        for count, sparql_phrase in enumerate(sparql_phrase_list):
            updated_sparql_phrase_list.append(f'({sparql_phrase})')
            if count < len(sparql_phrase_list)-1:
                updated_sparql_phrase_list.append(combinator_dict[combinator_list[count]])
            
        full_filter = " ".join(updated_sparql_phrase_list)
        return f'FILTER ({full_filter})'
        
    def create_query(self, source, destinations, values = None):
        """
        values is a dict with keys value, operation
        """
        if not isinstance(destinations, list):
            destinations = [destinations]
            
        #start prefix of quer
        query = []
        for key, val in self.namespaces.items():
            query.append(f'PREFIX {key}: <{val}>')
        for key, val in self.extra_namespaces.items():
            query.append(f'PREFIX {key}: <{val}>')
        
        #now for each destination, start adding the paths in the query
        all_triplets = {}
        for destination in destinations:
            triplets = self.get_shortest_path(source, destination, triples=True)
            all_triplets[destination] = triplets
        
        select_destinations = [f'?{self.strip_name(destination)}' for destination in destinations]
        query.append(f'SELECT DISTINCT {" ".join(select_destinations)}')
        query.append("WHERE {")
        
        #now add corresponding triples
        for destination in destinations:
            for triple in all_triplets[destination]:
                query.append("    ?%s %s ?%s ."%(self.strip_name(triple[0]), 
                                                 triple[1], 
                                                 self.strip_name(triple[2])))
        
        #now we have to add filters
        #filters are only needed if it is a dataproperty
        filter_text = ''
        
        if values is not None:
            lit_nodes = [node for node in self.g.nodes if 'node_type' in self.g.nodes[node].keys() and self.g.nodes[node]['node_type'] == 'literal']
            data_destinations = [destination for destination in destinations if destination in lit_nodes]
            if not len(data_destinations) == len(values):
                warnings.warn(f'Length of destinations and values are not same, found {len(data_destinations)} and {len(values)}')
                considered = " ".join(data_destinations[:len(values)])
                warnings.warn(f'Conditions are considered for {considered}')
            if len(values) > 0:
                filter_text = self.validate_values(data_destinations[:len(values)], values)
        query.append(filter_text)
        query.append('}')
        return '\n'.join(query)

        
            
            
            

        
        