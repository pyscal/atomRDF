from onto_network.term import OntoTerm 
from owlready2 import get_ontology

import os
import copy
import numpy as np
import itertools

class OntoParser:
    def __init__(self, infile, delimiter='/'):
        if os.path.exists(infile):
            self.tree = get_ontology(f'file://{infile}').load()
        else:
            raise FileNotFoundError(f'file {infile} not found!')
        self.attributes = {}
        self.attributes['class'] = {}
        self.attributes['object_property'] = {}
        self.attributes['data_property'] = {}        
        self.delimiter = delimiter
        self.classes = None
        self.namespaces = {self.tree.name: self.tree.base_iri}
        self.extra_namespaces = {'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        }
        self._parse_class()
        self._parse_object_property()
        self._parse_data_property()

    def __add__(self, ontoparser):
        """
        Add method; in principle it should add-
        - classes
        - attributes dict
        """
        for mainkey in ['class', 'object_property', 'data_property']:
            if mainkey in ontoparser.attributes.keys():
                for key, val in ontoparser.attributes[mainkey].items():
                    self.attributes[mainkey][key] = val        
        
        #now change classes
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

    def _strip_name(self, uri):
        uri_split = uri.split(self.delimiter)
        if len(uri_split)>1:
            return ":".join(uri_split[-2:])
        else:
            return uri
    
    def _strip_datatype(self, uri, delimiter='#'):
        uri_split = uri.split(delimiter)
        return uri_split[-1]        
        
    def _dict_to_lst(self, d):
        return [val for key, val in d.items()]
    
    def _get_subclasses(self, name):
        arg = self._in_which_bin_is_class(name)
        if arg is not None:
            return self.classes[arg]
        else:
            return [name]
                
    def _parse_data_property(self):
        for c in self.tree.data_properties():
            iri = c.iri
            dm = c.domain
            try:
                dm = [self._strip_name(d.iri) for d in dm[0].Classes]
            except:
                dm = [self._strip_name(d.iri) for d in dm]
            
            #now get subclasses
            dm = [self._get_subclasses(d) for d in dm]
            dm = list(itertools.chain(*dm))
            
            rn = c.range
            try:
                rn = [r.__name__ for r in rn[0].Classes if r is not None]
            except:
                rn = [r.__name__ for r in rn if r is not None]
            term = OntoTerm(iri)
            term.domain = dm
            term.range = rn
            term.node_type = 'data_property'
            self.attributes['data_property'][term.name] = term
            #assign this data
            for d in dm:
                if d!='07:owl#Thing':
                    self.attributes['class'][d].is_range_of.append(term.name)
            
    def _parse_object_property(self):
        for c in self.tree.object_properties():
            iri = c.iri
            dm = c.domain
            try:
                dm = [self._strip_name(d.iri) for d in dm[0].Classes]
            except:
                dm = [self._strip_name(d.iri) for d in dm]
            
            #now get subclasses
            dm = [self._get_subclasses(d) for d in dm]
            dm = list(itertools.chain(*dm))

            rn = c.range
            try:
                rn = [self._strip_name(r.iri) for r in rn[0].Classes]
            except:
                rn = [self._strip_name(r.iri) for r in rn]
            
            #now get subclasses
            rn = [self._get_subclasses(d) for d in rn]
            rn = list(itertools.chain(*rn))

            term = OntoTerm(iri)
            term.domain = dm
            term.range = rn
            term.node_type = 'object_property'
            self.attributes['object_property'][term.name] = term
            for d in dm:
                if d!='07:owl#Thing':
                    self.attributes['class'][d].is_range_of.append(term.name)
            for r in rn:
                if r!='07:owl#Thing':
                    self.attributes['class'][d].is_domain_of.append(term.name)
                
    def _parse_class_basic(self):
        classes = []
        for c in self.tree.classes():
            iri = c.iri
            #print(iri)
            try:
                subclasses = self.tree.search(subclass_of=getattr(self.tree, c.name))
                for sb in subclasses:
                    term = OntoTerm(sb.iri)
                    term.node_type ='class'
                    self.attributes['class'][term.name] = term
                subclasses = [self._strip_name(sb.iri) for sb in subclasses]
                classes.append(subclasses)
            except:
                term = OntoTerm(c.iri)
                term.node_type ='class'
                self.attributes['class'][term.name] = term                
                classes.append([self._strip_name(c.iri)])
        return classes
    
    def _aggregate_keys(self, dd):
        lst = copy.deepcopy(dd)
        #choose the first list
        large_list = []
        start = lst[0]
        #delete it from the main list
        nruns = len(lst)
        del lst[0]
        #now loop, if there is intersection add to this list
        while True:
            found = False
            index_to_delete = []
            for count, ls in enumerate(lst):
                common = len(list(set(start) & set(ls)))
                #print(common)
                if  common>0:
                    #common elements found! merge them
                    for l in ls:
                        start.append(l)
                    found = True
                    index_to_delete.append(count)
            if found:
                for ii in index_to_delete[::-1]:
                    del lst[ii]
            else:
                large_list.append(np.unique(start))
                if len(lst)==0:
                    break
                else:
                    start = lst[0]
                    del lst[0]
        return large_list
            
    def _parse_class(self):
        sub_classes = self._parse_class_basic()
        #now we have to go through and clean up sub classes
        sub_classes = self._aggregate_keys(sub_classes)
        self.classes = sub_classes
        
    def _in_which_bin_is_class(self, name):
        for count, lst in enumerate(self.classes):
            if name in lst:
                return count
        else:
            return None
    
            
            
        