"""
Graph module contains the basic RDFGraph object in atomrdf. This object gets a structure
as an input and annotates it with the CMSO ontology (PLDO and PODO too as needed). The annotated
object is stored in triplets.
"""

from rdflib import Graph, Literal,  XSD, RDF, RDFS, BNode, URIRef

import os
import numpy as np
import inspect
from ase.io import write
import copy
import pandas as pd
import yaml
import uuid
import json
import shutil
import tarfile

#from pyscal3.core import System
from pyscal3.atoms import Atoms

from atomrdf.visualize import visualize_graph
from atomrdf.network.network import OntologyNetwork
from atomrdf.network.ontology import read_ontology
from atomrdf.structure import System
import atomrdf.properties as prp
from atomrdf.stores import create_store
import atomrdf.json_io as json_io

from atomrdf.namespace import Namespace, CMSO, PLDO, PODO, ASMO

#read element data file
file_location = os.path.dirname(__file__).split('/')
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__),  'data/element.yml')
with open(file_location, 'r') as fin:
    element_indetifiers = yaml.safe_load(fin)


defstyledict = {
    "edgecolor": "#263238",
    "BNode": {"color": "#D9D9D9", 
              "shape": "box", 
              "style": "filled",
              "fontsize": "8",
              "fontname": "Helvetica"},
    "URIRef": {"color": "#C9DAF8", 
               "shape": "box", 
               "style": "filled",
               "fontsize": "8",
               "fontname": "Helvetica"},
    "Literal": {"color": "#E6B8AF", 
                "shape": "parallelogram", 
                "style": "filled",
                "fontsize": "8",
                "fontname": "Helvetica"},
}


def _replace_keys(refdict, indict):
    for key, val in indict.items():
        if key in refdict.keys():
            if isinstance(val, dict):
                _replace_keys(refdict[key], indict[key])
            else:
                refdict[key] = val
    return refdict


class KnowledgeGraph:
    def __init__(self, graph_file=None, 
        store="Memory", 
        store_file=None,
        identifier="http://default_graph",
        ontology=None,
        structure_store=None):
        

        create_store(self, store, identifier, 
            store_file=store_file,
            structure_store=structure_store)
        
        #start the storage

        #start binding
        self.graph.bind("cmso", CMSO)
        self.graph.bind("pldo", PLDO)
        
        if graph_file is not None:
            if os.path.exists(graph_file):
                self.graph.parse(graph_file)
        
        self.sample = None
        self.material = None
        self.sysdict = None
        self.sgraph = None
        if ontology is None:
            ontology = read_ontology()
        self.ontology = ontology
        self.terms = self.ontology.terms
        self._atom_ids = None
        self.store = store

    
    def add_structure(self, structure):
        structure.graph = self
        structure.to_graph()

    def _is_ontoterm(self, term):
        return type(term).__name__ == 'OntoTerm'

    def _modify_triple(self, triple):
        modified_triple = []
        for term in triple:
            if self._is_ontoterm(term):
                modified_triple.append(term.namespace_object)
            else:
                modified_triple.append(term)
        return tuple(modified_triple)

    def _check_domain_if_uriref(self, triple):
        found = True
        dm = self.value(triple[0], RDF.type)
        if dm is not None:
            #we need to check
            domain = triple[1].domain                
            if len(domain) > 0:
                if 'owl:Thing' not in domain:
                    found = False
                    for d in domain:
                        if d.split(':')[-1] in dm:
                            found = True
                            break
        return found
    
    def _check_domain(self, triple):
        if self._is_ontoterm(triple[1]):
            #check if type was provided
            found = True            
            if type(triple[0]).__name__ == 'URIRef': 
                found = self._check_domain_if_uriref(triple)

            if not found:
                raise ValueError(f'{dm} not in domain {domain} of {triple[1].name}')
            
            print(f'checked {triple[1].name}/{dm}')

    def add(self, triple, validate=True):
        """
        Force assumes that you are passing rdflib terms, defined with
        RDFLib Namespace
        """
        modified_triple = self._modify_triple(triple)

        #now we should validate if needed
        #check domain
        if validate:
            self._check_domain(triple)

        if str(modified_triple[2].toPython()) != 'None':
            self.graph.add(modified_triple)

    
    def triples(self, triple):
        modified_triple = self._modify_triple(triple)
        return self.graph.triples(modified_triple)

    def value(self, arg1, arg2):
        modified_double = self._modify_triple((arg1, arg2))
        return self.graph.value(modified_double[0], modified_double[1])

    def remove(self, triple):
        modified_triple = self._modify_triple(triple)
        return self.graph.remove(modified_triple)


    def _initialize_graph(self):
        """
        Create the RDF Graph from the data stored

        Parameters
        ----------
        names: bool
            if True, alphanumeric names will be used instead of random BNodes

        name_index: string
            Prefix to be added to identifiers, default 01        
        
        Returns
        -------
        None
        """        
        #extra triples
        self.add((CMSO.SimulationCellLength, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.LatticeParameter, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.Length, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM")))
        
        self.add((CMSO.SimulationCellAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.LatticeAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.Angle, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/DEG")))
        
        self.add((CMSO.LatticeVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.SimulationCellVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.PositionVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.Vector, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM")))
        

    def add_calculated_quantity(self, sample, propertyname, value, unit=None):
        prop = URIRef(f'{self._name}_{propertyname}')
        self.add((sample, CMSO.hasCalculatedProperty, prop))
        self.add((prop, RDF.type, CMSO.CalculatedProperty))
        self.add((prop, RDFS.label, Literal(propertyname)))
        self.add((prop, CMSO.hasValue, Literal(value)))
        if unit is not None:
            self.add((prop, CMSO.hasUnit, URIRef(f'http://qudt.org/vocab/unit/{unit}')))


    def inspect_sample(self):
        natoms = self.value(sample, CMSO.hasNumberOfAtoms).toPython()
        material = list([k[2] for k in self.triples((sample, CMSO.hasMaterial, None))])[0]
        defects = list([k[2] for k in self.triples((material, CMSO.hasDefect, None))])
        composition = list([k[2].toPython() for k in self.triples((material, CMSO.hasElementRatio, None))])
        crystalstructure = self.value(material, CMSO.hasStructure)
        spacegroupsymbol = self.value(crystalstructure, CMSO.hasSpaceGroupSymbol).toPython()

        lattice = self.value(sample, CMSO.hasNumberOfAtoms).toPython()
        defect_types = list([self.value(d, RDF.type).toPython() for d in defects])
        prop_nodes = list([k[2] for k in self.triples((sample, CMSO.hasCalculatedProperty, None))])
        props = list([self.value(prop_node, RDFS.label) for prop_node in prop_nodes])
        propvals = list([self.value(d, CMSO.hasValue).toPython() for d in prop_nodes])
        units = list([self.value(d, CMSO.hasUnit).toPython() for d in prop_nodes])
        st = []
        st.append(f'Sample with {natoms} atoms.\n')
        st.append("Material:\n")
        st.append(" ".join(composition))
        st.append("\n")
        st.append(f'Space Group symbol: {spacegroupsymbol}\n')
        if len(defect_types) > 0:
            st.append('With defects:\n')
            for d in defect_types:
                st.append(f'{d}\n')
        if len(props) > 0:
            st.append('With calculated properties:\n')
            for x in range(len(props)):
                st.append(f'{props[x]} with value: {propvals[x]} and unit: {units[x]}\n')

        return " ".join(st)

    def visualize(self, *args, **kwargs):
        """
        Vosualise the RDF tree of the Graph

        Parameters
        ----------
        backend: string, {'ipycytoscape', 'graphviz'}
            Chooses the backend with which the graph will be plotted. ipycytoscape provides an interactive, 
            but slow visualisation, whereas graphviz provides a non-interactive fast visualisation.

        edge_color: string
            Edge color of the boxes

        styledict: dict
            If provided, allows customisation of color and other properties.

        graph_attr: dict
            further attributes that allow customisation of graphs

        layoutname: string
            name of the layout for graph

        Returns
        -------

        Notes
        -----
        styledict has the following options. Refer to graphviz and ipycytoscape
        documentation for more details
        BNode:  
          color:  
          shape:  
          style:  
        URIRef:  
          color:  
          shape:  
          style:  
        Literal:  
          color:  
          shape:  
          style:         
        """
        self.visualise(*args, **kwargs)
        
    def visualise(self,
                  styledict=None, 
                  rankdir='BT',
                  hide_types=False,
                  workflow_view=False,
                  size=None,
                  layout='neato'):
        """
        Vosualise the RDF tree of the Graph

        Parameters
        ----------
        edge_color: string
            Edge color of the boxes

        styledict: dict
            If provided, allows customisation of color and other properties.

        graph_attr: dict
            further attributes that allow customisation of graphs

        layoutname: string
            name of the layout for graph

        Returns
        -------

        Notes
        -----
        styledict has the following options. Refer to graphviz and ipycytoscape
        documentation for more details
        BNode:  
          color:  
          shape:  
          style:  
        URIRef:  
          color:  
          shape:  
          style:  
        Literal:  
          color:  
          shape:  
          style:         
        """
        if size is not None:
            size = f"{size[0]},{size[1]}"

        sdict = defstyledict.copy()
        if styledict is not None:
            sdict = _replace_keys(sdict, styledict)
        
        return visualize_graph(self.graph, 
                               styledict=sdict, 
                               rankdir=rankdir,
                               hide_types=hide_types,
                               workflow_view=workflow_view,
                               size=size,
                               layout=layout)
    
    
    def write(self, filename, format="json-ld"):
        """
        Write the serialised version of the graph to a file

        Parameters
        ----------
        filename: string
            name of output file

        format: string, {'turtle', 'xml', 'json-ld', 'ntriples', 'n3'}
            output format to be written to 

        Returns
        -------
        None
        """

        with open(filename, "w") as fout:
            fout.write(self.graph.serialize(format=format))
    
    def archive(self, package_name, format='turtle', compress=True):
        """
        Publish a dataset from graph including per atom quantities
        """
        #first step make a folder
        if os.path.exists(package_name):
            raise ValueError(f'{package_name} already exists')
        if compress:
            if os.path.exists(f'{package_name}.tar.gz'):
                raise ValueError(f'{package_name} tarball already exists')
        
        os.mkdir(package_name)
        structure_store = f'{package_name}/rdf_structure_store' 
        os.mkdir(structure_store)

        #now go through each sample, and copy the file, at the same time fix the paths
        for sample in self.samples:
            filepath = self.value(URIRef(f'{sample}_Position'), CMSO.hasPath).toPython()
            shutil.copy(filepath, structure_store)
            
            #now we have to remove the old path, and fix new
            for val in ['Position', 'Species']:
                self.remove((URIRef(f'{sample}_{val}'), CMSO.hasPath, None))
            
                #assign corrected path
                new_relpath = "/".join(['rdf_structure_store', filepath.split('/')[-1]])
                self.graph.add((URIRef(f'{sample}_{val}'), CMSO.hasPath, Literal(new_relpath, datatype=XSD.string)))

        triple_file = os.path.join(package_name, 'triples')
        self.write(triple_file, format=format)

        if compress:
            with tarfile.open(f'{package_name}.tar.gz', "w:gz") as tar:
                tar.add(package_name, arcname=os.path.basename(package_name))
            shutil.rmtree(package_name)


    @classmethod
    def unarchive(cls, package_name, compress=True, 
        store="Memory", 
        store_file=None,
        identifier="http://default_graph",
        ontology=None):
        if compress:
            package_base_name = ".".join(package_name.split(".")[:-2])
            with tarfile.open(package_name) as fin: 
                fin.extractall(".")
            #os.remove(package_name)
            #copy things out

        return cls(store=store, store_file=store_file,
            identifier=identifier, 
            graph_file=f'{package_base_name}/triples', 
            structure_store=f'{package_base_name}/rdf_structure_store',
            ontology=ontology)
    
    def query(self, inquery):
        """
        Query the graph using SPARQL

        Parameters
        ----------
        inquery: string
            SPARQL query to be executed

        Returns
        -------
        res: pandas DataFrame
            pandas dataframe results
        """
        res = self.graph.query(inquery)
        if res is not None:
            for line in inquery.split('\n'):
                if 'SELECT DISTINCT' in line:
                    break
            labels = [x[1:] for x in line.split()[2:]]
            return pd.DataFrame(res, columns=labels)
        raise ValueError("SPARQL query returned None")

    def auto_query(self, source, destination, 
        condition=None, 
        return_query=False, 
        enforce_types=None):

        if enforce_types is None:
            for val in [True, False]:
                query = self.ontology.create_query(source, destination, 
                    condition=condition, enforce_types=val)
                if return_query:
                    return query
                res = self.query(query)
                if len(res) != 0:
                    return res
        else:
            query = self.ontology.create_query(source, destination, 
                condition=condition, enforce_types=enforce_types)
            if return_query:
                return query
            res = self.query(query)

        return res    

    #################################
    # Methods to interact with sample
    #################################
    def query_sample(self, destination, condition=None, return_query=False, enforce_types=None):
        return self.auto_query(self.ontology.terms.cmso.AtomicScaleSample, destination,
            condition=condition, return_query=return_query, enforce_types=enforce_types)

    @property
    def n_samples(self):
        """
        Number of samples in the Graph
        """

        return len([x for x in self.triples((None, RDF.type, CMSO.AtomicScaleSample))])
    
    @property
    def samples(self):
        """
        Returns a list of all Samples in the graph
        """

        return [x[0] for x in self.triples((None, RDF.type, CMSO.AtomicScaleSample))]
        
    def iterate_graph(self, item, create_new_graph=False):
        if create_new_graph:
            self.sgraph = KnowledgeGraph()
        triples = list(self.triples((item, None, None)))
        for triple in triples:
            self.sgraph.graph.add(triple)
            self.iterate_graph(triple[2])
    
    def get_sample(self, sample, no_atoms=False):
        """
        Get the Sample as an RDFGraph

        Parameters
        ----------
        sample: string
            sample id

        no_atoms: bool, optional
            if True, returns the number of atoms in the sample

        Returns
        -------
        sgraph: :py:class:`RDFGraph`
            the RDFGraph of the queried sample

        """

        self.iterate_graph(sample, create_new_graph=True)
        if no_atoms:
            na = self.sgraph.value(sample, CMSO.hasNumberOfAtoms).toPython()
            return self.sgraph, na
        return self.sgraph
        
    def get_system_from_sample(self, sample):
        """
        Get a pyscal :py:class:`pyscal.core.System` from the selected sample

        Parameters
        ----------
        sample: string
            sample id

        Returns
        -------
        system: :py:class:`pyscal.core.System`
            corresponding system
        """

        simcell = self.value(sample, CMSO.hasSimulationCell)
        cell_vectors = [[], [], []]

        for s in self.triples((simcell, CMSO.hasVector, None)):
            cell_vectors[0].append(self.value(s[2], CMSO.hasComponent_x).toPython())
            cell_vectors[1].append(self.value(s[2], CMSO.hasComponent_y).toPython())
            cell_vectors[2].append(self.value(s[2], CMSO.hasComponent_z).toPython())
        
        #cell_vectors
        filepath = self.value(URIRef(f'{sample}_Position'), CMSO.hasPath).toPython()
        position_identifier = self.value(URIRef(f'{sample}_Position'), CMSO.hasIdentifier).toPython()
        species_identifier = self.value(URIRef(f'{sample}_Species'), CMSO.hasIdentifier).toPython()

        #open the file for reading
        with open(filepath, 'r') as fin:
            data = json.load(fin)
            positions = data[position_identifier]['value']
            species = data[species_identifier]['value']

        atoms = {"positions": positions, "species": species}
        at = Atoms()
        at.from_dict(atoms)
        sys = System()
        sys.box = cell_vectors
        sys.atoms = at       
        return sys

    def to_file(self, sample, filename=None, format="lammps-dump"):
        """
        Save a given sample to a file

        Parameters
        ----------
        sample
            ID of the sample

        filename: string
            name of output file

        format: string, {"lammps-dump","lammps-data", "poscar"}

        Returns
        -------
        None
        """

        if filename is None:
            filename = os.path.join(os.getcwd(), "out")
        
        sys = self.get_system_from_sample(sample)
        
        if format=="ase":
            return sys.write.ase()
        elif format=='poscar':
            asesys = sys.write.ase()
            write(filename, asesys, format="vasp")
        else:
            asesys = sys.write.ase()
            write(filename, asesys, format=format)
