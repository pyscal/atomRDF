import networkx as nx
import graphviz
import matplotlib.pyplot as plt
import numpy as np
import os
from pyscal_rdf.parser import OntoParser

owlfile = os.path.join(os.path.dirname(__file__), "data/cmso.owl")

def _replace_name(name):
    return ".".join(name.split(':'))

class OntologyNetwork:
    """
    Network representation of Onto
    """
    def __init__(self, infile=None):
        if infile is None:
            infile = owlfile
            
        self.g = nx.DiGraph()
        self.onto = OntoParser(infile)
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

    @property
    def attributes(self):
        return self.onto.attributes

    def __radd__(self, ontonetwork):
        return self.__add__(ontonetwork)

    def get_shortest_path(self, source, target):
        path = nx.shortest_path(self.g, source=source, target=target)
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
        path = self.get_shortest_path(source="Sample", target=target)
        triplets = []
        for x in range(len(path)//2):
            triplets.append(path[2*x:2*x+3])
        return triplets
        
    def formulate_query(self, target, value):
        #first get triplets
        triplets = self.get_path_from_sample(target)
        #start building query
        query = self._formulate_query_path(triplets)
        query.append(self._formulate_filter_expression(triplets, value))
        query.append("}")
        query = " ".join(query)
        return query
        
    
    def _formulate_query_path(self, triplets):
        query = []
        query.append("PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>")
        query.append("PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>")
        query.append("PREFIX podo: <https://purls.helmholtz-metadaten.de/podo/>")
        query.append("PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>")
        query.append("SELECT DISTINCT ?sample")
        query.append("WHERE {")
        for triple in triplets:
            query.append("    ?%s %s ?%s ."%(triple[0].lower(), 
                                                  triple[1], 
                                                  triple[2].lower()))
        return query
    
    def _formulate_filter_expression(self, triplets, value):                       
        value, datatype = self._check_value(value)      
        last_val = self.g.nodes[triplets[-1][-1]]
        last_val_name = triplets[-1][-1].lower()
        
        #if it is nodetype data
        if last_val['node_type'] == "data":
            if datatype == "multi_string":
                qstr = self._formulate_or_string_query(last_val, 
                                                   last_val_name, 
                                                   value)
            elif datatype == "multi_number":
                qstr = self._formulate_range_number_query(last_val, 
                                                   last_val_name, 
                                                   value)
            else:
                qstr = self._formulate_equal_query(last_val, 
                                                   last_val_name, 
                                                   value)
            return qstr
        else:
            raise NotImplementedError("Non-data queries are not implemented")
    
    def _check_value(self, value):
        if isinstance(value, list):
            if not len(value) == 2:
                raise ValueError("value can be maximum length 2")
        else:
            value = [value]
        if all(isinstance(x, str) for x in value):
            datatype = "string"
        elif all(isinstance(x, (int, float)) for x in value):
            datatype = "number"
        else:
            raise TypeError("Values have to be of same type")
        if len(value) == 1:
            datatype = f'single_{datatype}'
        else:
            datatype = f'multi_{datatype}'
        return value, datatype
    
    
    def _formulate_equal_query(self, last_val, last_val_name, value):
        qstr = "FILTER (?%s=\"%s\"^^xsd:%s)"%(last_val_name, 
                                              str(value[0]), 
                                              last_val['dtype'])
        return qstr
    
    def _formulate_or_string_query(self, last_val, last_val_name, value):
        qstr = "FILTER (?%s=\"%s\"^^xsd:%s || ?%s=\"%s\"^^xsd:%s)"%(last_val_name, 
                                                                    str(value[0]), 
                                                                    last_val['dtype'],
                                                                    last_val_name, 
                                                                    str(value[1]), 
                                                                    last_val['dtype'],)
        return qstr
    
    def _formulate_range_number_query(self, last_val, last_val_name, value):
        value = np.sort(value)
        qstr = "FILTER (?%s >= \"%s\"^^xsd:%s && ?%s <= \"%s\"^^xsd:%s)"%(last_val_name, 
                                                                    str(value[0]), 
                                                                    last_val['dtype'],
                                                                    last_val_name, 
                                                                    str(value[1]), 
                                                                    last_val['dtype'],)
        return qstr

        
            
            
            

        
        