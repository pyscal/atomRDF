import os
import copy
import numpy as np
import itertools

from atomrdf.network.term import OntoTerm, strip_name
from atomrdf.network.patch import patch_terms
from rdflib import Graph, RDF, RDFS, OWL, BNode, URIRef

class OntoParser:
    def __init__(self, infile, format='xml'):
        if not os.path.exists(infile):
            raise FileNotFoundError(f"file {infile} not found!")
        
        self.graph = Graph()
        self.graph.parse(infile, format=format)
        self.classes = []
        self.attributes = {}
        self.attributes["class"] = {}
        self.attributes["object_property"] = {}
        self.attributes["data_property"] = {}
        self.attributes["data_nodes"] = {}
        self.mappings = {}
        self.namespaces = {}
        self.extra_namespaces = {}

        self.extract_classes()
        self.extract_relations(relation_type="union")
        self.extract_relations(relation_type="intersection")
        self.add_classes_to_attributes()
        self.parse_subclasses()
        self.parse_equivalents()
        self.parse_named_individuals()
        self.extract_object_properties()
        self.extract_data_properties()
        self.recheck_namespaces()        

    def __add__(self, ontoparser):
        """
        Add method; in principle it should add-
        - classes
        - attributes dict
        """
        for mainkey in ["class", "object_property", "data_property"]:
            if mainkey in ontoparser.attributes.keys():
                for key, val in ontoparser.attributes[mainkey].items():
                    self.attributes[mainkey][key] = val

        # now change classes
        if ontoparser.classes is not None:
            for clx in ontoparser.classes:
                self.classes.append(clx)

        for key, val in ontoparser.namespaces.items():
            self.namespaces[key] = val

        for key, val in ontoparser.extra_namespaces.items():
            self.extra_namespaces[key] = val

        return self

    def __radd__(self, ontoparser):
        return self.__add__(ontoparser)
    
    @property
    def base_iri(self):
        base_iri = None
        for s in self.graph.subjects(RDF.type, OWL.Ontology):
            base_iri = str(s)
        return base_iri
    


    def recheck_namespaces(self):
        for mainkey in ["class", "object_property", "data_property"]:
            for key, val in self.attributes[mainkey].items():
                namespace = self.attributes[mainkey][key].namespace
                if namespace not in self.namespaces.keys():
                    self.namespaces[namespace] = self.attributes[mainkey][
                        key
                    ].namespace_with_prefix

    def extract_classes(self):
        self.classes = list(self.graph.subjects(RDF.type, OWL.Class))

    def extract_object_properties(self):
        object_properties = list(self.graph.subjects(RDF.type, OWL.ObjectProperty))
        for cls in object_properties:
            term = self.create_term(cls)
            term.domain = self.get_domain(cls)
            term.range = self.get_range(cls)
            term.node_type = "object_property"
            self.attributes["object_property"][term.name] = term            

    def extract_data_properties(self):
        data_properties = list(self.graph.subjects(RDF.type, OWL.DatatypeProperty))
        for cls in data_properties:
            term = self.create_term(cls)
            term.domain = self.get_domain(cls)
            rrange = self.get_range(cls)
            rrange = [x.split(":")[-1] for x in rrange]
            rrange = patch_terms(term.uri, rrange)

            term.range = rrange
            term.node_type = "data_property"
            self.attributes["data_property"][term.name] = term

            #now create data nodes
            data_term = OntoTerm()
            data_term.name = term.name + "value"
            data_term.node_type = "data_node"
            self.attributes["data_property"][term.name].associated_data_node = data_term.name
            self.attributes["data_nodes"][data_term.name] = data_term   
    
    def extract_values(self, subject, predicate):
        vallist = list([x[2] for x in self.graph.triples((subject, predicate, None))])
        if len(vallist) > 0:
            return vallist[0]
        else:
            return None
            
    def extract_relations(self, relation_type):
        if relation_type == "union":
            owl_term = OWL.unionOf
        elif relation_type == "intersection":
            owl_term = OWL.intersectionOf
            
        to_delete = []
        for term in self.classes:
            if isinstance(term, BNode):
                union_term = self.extract_values(term, owl_term)
                if union_term is not None:
                    unravel_list = []
                    self.unravel_relation(union_term, unravel_list)
                    self.mappings[term.toPython()] = {"type": relation_type, "items": [strip_name(item.toPython()) for item in unravel_list]}
                    to_delete.append(term)
        for term in to_delete:
            self.classes.remove(term)

    def add_classes_to_attributes(self):
        for cls in self.classes:
            term = self.create_term(cls)
            term.node_type = "class"
            self.attributes["class"][term.name] = term

    def get_description(self, cls):
        comment = self.graph.value(cls, URIRef('http://purl.obolibrary.org/obo/IAO_0000115'))
        if comment is None:
            comment = self.graph.value(cls, URIRef('http://www.w3.org/2000/01/rdf-schema#comment'))
        if comment is None:
            comment = ""
        return comment

    def lookup_node(self, term):
        if isinstance(term, BNode):
            #lookup needed
            term_name = term.toPython()
            if term_name in self.mappings:
                terms = self.mappings[term_name]['items']
            else:
                terms = [strip_name(term.toPython())]
        else:
            terms = [strip_name(term.toPython())]
        #so here we map the domain and range wrt to other heirarchies
        additional_terms = []
        #first get subclasses which will share the domain and range
        for term in terms:
            #check if such a thing exists in the class 
            if term in self.attributes['class']:
                #get the subclasses
                additional_terms += self.attributes['class'][term].subclasses
                #get the equivalent classes
                additional_terms += self.attributes['class'][term].equivalent_classes
                #get the named individuals
                additional_terms += self.attributes['class'][term].named_individuals
        #add additiona terms to terms
        terms += additional_terms
        return terms

    def lookup_class(self, term):
        if isinstance(term, BNode):
            term = term.toPython()
        else:
            term = strip_name(term.toPython())
        #print(term)
        if term in self.attributes['class']:
            return [self.attributes['class'][term].name]
        elif term in self.mappings:
            return self.mappings[term]['items']

    def get_domain(self, cls):
        domain = []
        for triple in self.graph.triples((cls, URIRef('http://www.w3.org/2000/01/rdf-schema#domain'), None)):
            domain_term = self.lookup_node(triple[2])
            for term in domain_term:
                domain.append(term)
        return domain

    def get_range(self, cls):
        rrange = []
        for triple in self.graph.triples((cls, URIRef('http://www.w3.org/2000/01/rdf-schema#range'), None)):
            range_term = self.lookup_node(triple[2])
            for term in range_term:
                rrange.append(term)
        return rrange
        
    def create_term(self, cls):
        iri = cls.toPython()
        term = OntoTerm(iri)
        term.description = self.get_description(cls)
        term._object = cls
        return term
                    
    def unravel_relation(self, term, unravel_list):
        if term == RDF.nil:
            return 
        first_term = self.graph.value(term, RDF.first)
        if first_term not in unravel_list:
            unravel_list.append(first_term)
        second_term =  self.graph.value(term, RDF.rest)
        self.unravel_relation(second_term, unravel_list)

    def parse_subclasses(self):
        for key, cls in self.attributes['class'].items():
            for triple in self.graph.triples((cls._object, RDFS.subClassOf, None)):
                superclasses = self.lookup_class(triple[2])
                for superclass in superclasses:
                    self.attributes['class'][superclass].subclasses.append(cls.name)
    
    def parse_equivalents(self):
        for key, cls in self.attributes['class'].items():
            for triple in self.graph.triples((cls._object, OWL.equivalentClass, None)):
                equivalent = triple[2]
                self.attributes['class'][strip_name(equivalent)].equivalent_classes.append(cls.name)
                cls.equivalent_classes.append(strip_name(equivalent))
    
    def parse_named_individuals(self):
        named_individuals = list(self.graph.subjects(RDF.type, OWL.NamedIndividual))
        for cls in named_individuals:
            #find parent
            term = self.create_term(cls)
            self.attributes["class"][term.name] = term
            parents = list(self.graph.objects(cls, RDF.type))
            for parent in parents:
                if parent != OWL.NamedIndividual:
                    self.attributes["class"][strip_name(parent.toPython())].named_individuals.append(term.name)
