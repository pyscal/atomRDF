"""
Graph module contains the basic RDFGraph object in atomrdf. This object gets a structure
as an input and annotates it with the CMSO ontology (PLDO and PODO too as needed). The annotated
object is stored in triplets.

NOTES
-----
- To ensure domain and range checking works as expected, always add type before adding further properties!

Classes
-------
- KnowledgeGraph: Represents a knowledge graph that stores and annotates structure objects.

Attributes
----------
- defstyledict: A dictionary containing default styles for visualizing the graph.

"""


from rdflib import Graph, XSD, RDF, RDFS, BNode, URIRef

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
import logging
import warnings
import re
import pickle

# from pyscal3.core import System
from pyscal3.atoms import Atoms

from atomrdf.visualize import visualize_graph, visualize_provenance
from atomrdf.network.network import OntologyNetwork
from atomrdf.network.ontology import read_ontology
from atomrdf.structure import System
import atomrdf.properties as prp
from atomrdf.stores import create_store
import atomrdf.json_io as json_io
from atomrdf.workflow.workflow import Workflow
from atomrdf.sample import Sample
import atomrdf.mp as amp 

from atomrdf.namespace import Namespace, CMSO, PLDO, PODO, ASMO, PROV, MATH, UNSAFECMSO, UNSAFEASMO, Literal

# read element data file
file_location = os.path.dirname(__file__).split("/")
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__), "data/element.yml")
with open(file_location, "r") as fin:
    element_indetifiers = yaml.safe_load(fin)


defstyledict = {
    "edgecolor": "#263238",
    "BNode": {
        "color": "#D9D9D9",
        "shape": "box",
        "style": "filled",
        "fontsize": "8",
        "fontname": "Helvetica",
    },
    "URIRef": {
        "color": "#C9DAF8",
        "shape": "box",
        "style": "filled",
        "fontsize": "8",
        "fontname": "Helvetica",
    },
    "Literal": {
        "color": "#E6B8AF",
        "shape": "parallelogram",
        "style": "filled",
        "fontsize": "8",
        "fontname": "Helvetica",
    },
}

def _clean_string(input_string):
    input_string = re.sub(r'\W', '_', input_string)
    if input_string[0].isdigit():
        input_string = "s" + input_string
    return input_string

def _replace_keys(refdict, indict):
    for key, val in indict.items():
        if key in refdict.keys():
            if isinstance(val, dict):
                _replace_keys(refdict[key], indict[key])
            else:
                refdict[key] = val
    return refdict


def _dummy_log(str):
    pass

def _name(term):
    try:
        return str(term.toPython())
    except:
        return str(term)
    
def _prepare_log(file):
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(file)
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    return logger

class KnowledgeGraph:
    """
    Represents a knowledge graph.

    Parameters
    ----------
    graph_file : str, optional
        The path to the graph file to be parsed. Default is None.
    store : str, optional
        The type of store to use. Default is "Memory".
    store_file : str, optional
        The path to the store file. Default is None.
    identifier : str, optional
        The identifier for the graph. Default is "http://default_graph".
    ontology : Ontology, optional
        The ontology object to be used. Default is None.
    structure_store : StructureStore, optional
        The structure store object to be used. Default is None.
    enable_log : bool, optional
        Whether to enable logging. Default is False.
        If true, a log file named atomrdf.log will be created in the current working directory.

    Attributes
    ----------
    graph : rdflib.Graph
        The RDF graph.
    sgraph : rdflib.Graph
        The structure graph for a single chosen sample
    ontology : Ontology
        The ontology object.
    terms : dict
        The dictionary of ontology terms.
    store : str
        The type of store used.

    Methods
    -------
    add_structure(structure)
        Add a structure to the knowledge graph.
    add(triple, validate=True)
        Add a triple to the knowledge graph.
    triples(triple)
        Return the triples in the knowledge graph that match the given triple pattern.
    """

    def __init__(
        self,
        graph_file=None,
        store="Memory",
        store_file=None,
        identifier="http://default_graph",
        ontology=None,
        structure_store=None,
        enable_log=False,
    ):
        """
        Initialize the KnowledgeGraph object.

        Parameters
        ----------
        graph_file : str, optional
            The path to the graph file to be parsed. Default is None.
        store : str, optional
            The type of store to use. Default is "Memory".
        store_file : str, optional
            The path to the store file. Default is None.
        identifier : str, optional
            The identifier for the graph. Default is "http://default_graph".
        ontology : Ontology, optional
            The ontology object to be used. Default is None.
        structure_store : StructureStore, optional
            The structure store object to be used. Default is None.
        enable_log : bool, optional
            Whether to enable logging. Default is False.
        """

        create_store(
            self,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )

        # enable logging
        if enable_log:
            logger = _prepare_log(os.path.join(os.getcwd(), "atomrdf.log"))
            self.log = logger.info
        else:
            self.log = _dummy_log

        # start binding
        self.graph.bind("cmso", CMSO)
        self.graph.bind("pldo", PLDO)

        if graph_file is not None:
            if os.path.exists(graph_file):
                self.graph.parse(graph_file)

        self.sgraph = None
        if ontology is None:
            ontology = read_ontology()
        self.ontology = ontology
        self.terms = self.ontology.terms
        self.store = store
        self._n_triples = 0
        self._initialize_graph()
        self.workflow = Workflow(self)

    def add_structure(self, structure):
        """
        Add a structure to the knowledge graph.

        Parameters
        ----------
        structure : Structure
            The structure object to be added.

        Returns
        -------
        None

        Notes
        -----
        This method adds a structure object to the knowledge graph. The structure object should be an instance of the Structure class.
        The structure object is assigned to the graph and converted to RDF format.
        """
        structure.graph = self
        structure.to_graph()

    def _is_valid(self, input_list):
        valid = False
        flat_list = []
        for x in input_list:
            if isinstance(x, list):
                flat_list.extend(x)
            else:
                flat_list.append(x)
        for x in flat_list:
            if x is not None:
                valid = True
                break
        return valid

    def _is_ontoterm(self, term):
        return type(term).__name__ == "OntoTerm"

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
            # we need to check
            domain = triple[1].domain
            if len(domain) > 0:
                if "owl:Thing" not in domain:
                    if triple[1].namespace_with_prefix not in dm:
                        # cross ontology term
                        self.log(
                            f"ignoring possible cross ontology connection between {triple[1].namespace} and {dm}"
                        )
                        return True, None

                    found = False
                    for d in domain:
                        if d.split(":")[-1] in dm:
                            found = True
                            break
        return found, dm

    #    def _check_domain_if_ontoterm(self, triple):
    #        found = True
    #        domain = triple[0].domain
    #        if len(domain) > 0:
    #            if 'owl:Thing' not in domain:
    #                if triple[1].namespace != triple[0].namespace:
    #                    #cross ontology term
    #                    self.log(f'ignoring possible cross ontology connection between {triple[1].namespace} and {triple[0].namespace}')
    #                    return True, None
    #                found = False
    #                if triple[1].name in domain:
    #                    found = True
    #        return found, triple[0].name

    def _check_domain(self, triple):
        if self._is_ontoterm(triple[1]):
            # check if type was provided
            found = True
            dm = None

            if type(triple[0]).__name__ == "URIRef":
                found, dm = self._check_domain_if_uriref(triple)

            # elif self._is_ontoterm(triple[0]):
            #    found, dm = self._check_domain_if_ontoterm(triple)

            if not found:
                raise ValueError(f"{dm} not in domain of {triple[1].name}")

            self.log(f"checked {triple[1].name} against domain {dm}")

    def _check_range_if_uriref(self, triple):
        found = True
        rn = self.value(triple[2], RDF.type)

        if rn is not None:
            # we need to check
            rang = triple[1].range
            if len(rang) > 0:
                if "owl:Thing" not in rang:
                    if triple[1].namespace_with_prefix not in rn:
                        # cross ontology term
                        self.log(
                            f"ignoring possible cross ontology connection between {triple[1].namespace} and {rn}"
                        )
                        return True, None

                    found = False
                    for r in rang:
                        if r.split(":")[-1] in rn:
                            found = True
                            break
        return found, rn

    #    def _check_range_if_ontoterm(self, triple):
    #        found = True
    #        rang = triple[1].range
    #        if len(rang) > 0:
    #            if 'owl:Thing' not in rang:
    #                if triple[1].namespace != triple[2].namespace:
    #                    #cross ontology term
    #                    self.log(f'ignoring possible cross ontology connection between {triple[1].namespace} and {triple[2].namespace}')
    #                    return True, None
    #
    #                found = False
    #                if triple[2].name in rang:
    #                    found = True
    #        return found, triple[2].name

    def _check_range_if_literal(self, triple):
        found = True
        if triple[2].datatype is None:
            self.log(
                f"WARNING: {triple[1].name} has a range with unspecified datatype!"
            )
            warnings.warn(f"{triple[1].name} has a range with unspecified datatype!")
            return True, None

        destination_range = triple[2].datatype.toPython().split("#")[-1]

        if destination_range == "string":
            destination_range = "str"
        elif destination_range == "integer":
            destination_range = "int"

        rang = triple[1].range
        if len(rang) > 0:
            found = False
            if destination_range in rang:
                found = True
        return found, destination_range

    def _check_range(self, triple):
        if self._is_ontoterm(triple[1]):
            # check if type was provided
            found = True
            dm = None

            if type(triple[2]).__name__ == "URIRef":
                found, dm = self._check_range_if_uriref(triple)

            # elif self._is_ontoterm(triple[2]):
            #    found, dm = self._check_range_if_ontoterm(triple)

            elif type(triple[2]).__name__ == "Literal":
                found, dm = self._check_range_if_literal(triple)

            if not found:
                raise ValueError(f"{dm} not in range of {triple[1].name}")

            self.log(f"checked {triple[1].name} against range {dm}")

    def add(self, triple, validate=True):
        """
        Add a triple to the knowledge graph.

        Parameters
        ----------
        triple : tuple
            The triple to be added in the form (subject, predicate, object).
        validate : bool, optional
            Whether to validate the triple against the domain and range. Default is True.

        Returns
        -------
        None

        Notes
        -----
        This method adds a triple to the knowledge graph. The triple should be provided as a tuple in the form (subject, predicate, object).
        By default, the triple is validated against the domain and range. If the `validate` parameter is set to False, the validation is skipped.

        Examples
        --------
        >>> graph = Graph()
        >>> graph.add(("Alice", "likes", "Bob"))
        >>> graph.add(("Bob", "age", 25), validate=False)
        """
        modified_triple = self._modify_triple(triple)

        self.log(f"attempting to add triple: {self._n_triples}")
        self.log(f"- {modified_triple[0].toPython()}")
        self.log(f"- {modified_triple[1].toPython()}")
        self.log(f"- {modified_triple[2].toPython()}")

        if validate:
            self._check_domain(triple)
            self._check_range(triple)

        if str(modified_triple[2].toPython()) == "None":
            self.log(f"rejecting None valued triple")
            return

        self.graph.add(modified_triple)
        self._n_triples += 1

        self.log("added")

    def triples(self, triple):
        """
        Return the triples in the knowledge graph that match the given triple pattern.

        Parameters
        ----------
        triple : tuple
            The triple pattern to match in the form (subject, predicate, object).

        Returns
        -------
        generator
            A generator that yields the matching triples.

        Examples
        --------
        >>> graph = KnowledgeGraph()
        >>> graph.add(("Alice", "likes", "Bob"))
        >>> graph.add(("Alice", "dislikes", "Charlie"))
        >>> graph.add(("Bob", "likes", "Alice"))
        >>> for triple in graph.triples(("Alice", None, None)):
        ...     print(triple)
        ('Alice', 'likes', 'Bob')
        ('Alice', 'dislikes', 'Charlie')
        """
        modified_triple = self._modify_triple(triple)
        return self.graph.triples(modified_triple)

    def value(self, arg1, arg2):
        """
        Get the value of a triple in the knowledge graph.

        Parameters
        ----------
        arg1 : object
            The subject of the triple.
        arg2 : object
            The predicate of the triple.

        Returns
        -------
        object or None
            The value of the triple if it exists, otherwise None.

        Notes
        -----
        This method retrieves the value of a triple in the knowledge graph. The triple is specified by providing the subject and predicate as arguments. 
        If the triple exists in the graph, the corresponding value is returned. If the triple does not exist, None is returned.

        Examples
        --------
        >>> graph = KnowledgeGraph()
        >>> graph.add(("Alice", "likes", "Bob"))
        >>> value = graph.value("Alice", "likes")
        >>> print(value)
        Bob
        """
        modified_double = self._modify_triple((arg1, arg2))
        return self.graph.value(modified_double[0], modified_double[1])

    def remove(self, triple):
        """
        Remove a triple from the knowledge graph.

        Parameters
        ----------
        triple : tuple
            The triple to be removed in the form (subject, predicate, object).

        Returns
        -------
        None

        Notes
        -----
        This method removes a triple from the knowledge graph. The triple should be provided as a tuple in the form (subject, predicate, object).

        Examples
        --------
        >>> graph = KnowledgeGraph()
        >>> graph.add(("Alice", "likes", "Bob"))
        >>> graph.remove(("Alice", "likes", "Bob"))
        """
        modified_triple = self._modify_triple(triple)
        return self.graph.remove(modified_triple)

    def create_node(self, namestring, classtype, label=None):
        """
        Create a new node in the graph.

        Parameters
        ----------
        namestring : str
            The name of the node.
        classtype : Object from a given ontology
            The class type of the node.

        Returns
        -------
        URIRef
            The newly created node.

        """
        item = URIRef(namestring)
        self.add((item, RDF.type, classtype))
        if label is not None:
            self.add((item, RDFS.label, Literal(_clean_string(label))))
        return item

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
        # extra triples
        self.add((CMSO.SimulationCellLength, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.LatticeParameter, RDFS.subClassOf, CMSO.Length))
        self.add(
            (CMSO.Length, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM"))
        )

        self.add((CMSO.SimulationCellAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.LatticeAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.Angle, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/DEG")))

        self.add((CMSO.LatticeVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.SimulationCellVector, RDFS.subClassOf, CMSO.Vector))
        # self.add((CMSO.PositionVector, RDFS.subClassOf, CMSO.Vector))
        self.add(
            (CMSO.Vector, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM"))
        )

    def add_calculated_quantity(self, sample, propertyname, value, unit=None):
        """
        Add a calculated quantity to a sample.

        Parameters
        ----------
        sample : URIRef
            The URIRef of the sample to which the calculated quantity is being added.
        propertyname : str
            The name of the calculated property.
        value : str
            The value of the calculated property.
        unit : str, optional
            The unit of the calculated property. Default is None.
            The unit should be from QUDT. See http://qudt.org/vocab/unit/

        Returns
        -------
        None

        Notes
        -----
        This method adds a calculated quantity to a sample in the knowledge graph. The calculated quantity is represented as a triple with the sample as the subject, the calculated property as the predicate, and the value as the object. The calculated property is created as a node in the graph with the given name and value. If a unit is provided, it is also added as a property of the calculated property node.

        Examples
        --------
        >>> graph = KnowledgeGraph()
        >>> sample = graph.create_node("Sample1", CMSO.Sample)
        >>> graph.add_calculated_quantity(sample, "energy", "10.5", "eV")
        """
        prop = self.create_node(propertyname, CMSO.CalculatedProperty)
        self.add((sample, CMSO.hasCalculatedProperty, prop))
        self.add((prop, RDFS.label, Literal(propertyname)))
        self.add((prop, ASMO.hasValue, Literal(value)))
        if unit is not None:
            self.add((prop, ASMO.hasUnit, URIRef(f"http://qudt.org/vocab/unit/{unit}")))

    def inspect_sample(self, sample):
        """
        Inspects a sample and retrieves information about its atoms, material, defects, composition,
        crystal structure, space group, calculated properties, and units.

        Parameters
        ----------
        sample: The sample to inspect.

        Returns
        -------
        string: A string containing the information about the sample.

        """
        natoms = self.value(sample, CMSO.hasNumberOfAtoms).toPython()
        material = list([k[2] for k in self.triples((sample, CMSO.hasMaterial, None))])[
            0
        ]
        defects = list([k[2] for k in self.triples((material, CMSO.hasDefect, None))])
        composition = list(
            [
                k[2].toPython()
                for k in self.triples((material, CMSO.hasElementRatio, None))
            ]
        )
        crystalstructure = self.value(material, CMSO.hasStructure)
        spacegroupsymbol = self.value(crystalstructure, CMSO.hasSpaceGroupSymbol).toPython()

        lattice = self.value(sample, CMSO.hasNumberOfAtoms).toPython()
        defect_types = list([self.value(d, RDF.type).toPython() for d in defects])
        prop_nodes = list(
            [k[2] for k in self.triples((sample, CMSO.hasCalculatedProperty, None))]
        )
        props = list([self.value(prop_node, RDFS.label) for prop_node in prop_nodes])
        propvals = list([self.value(d, ASMO.hasValue).toPython() for d in prop_nodes])
        units = list([self.value(d, ASMO.hasUnit).toPython() for d in prop_nodes])
        st = []
        st.append(f"Sample with {natoms} atoms.\n")
        st.append("Material:\n")
        st.append(" ".join(composition))
        st.append("\n")
        st.append(f"Space Group symbol: {spacegroupsymbol}\n")
        if len(defect_types) > 0:
            st.append("With defects:\n")
            for d in defect_types:
                st.append(f"{d}\n")
        if len(props) > 0:
            st.append("With calculated properties:\n")
            for x in range(len(props)):
                st.append(
                    f"{props[x]} with value: {propvals[x]} and unit: {units[x]}\n"
                )

        return " ".join(st)

    def visualize(self, *args, **kwargs):
        """
        Visualizes the graph using the specified arguments.

        This method is a wrapper around the `visualise` method and passes the same arguments to it.

        Parameters
        ----------
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

        Returns
        -------
        dot: The visualization of the RDF tree.
        """
        self.visualise(*args, **kwargs)

    def visualise(
        self,
        styledict=None,
        rankdir="BT",
        hide_types=False,
        workflow_view=False,
        sample_view=False,
        size=None,
        layout="neato",
    ):
        """
        Visualize the RDF tree of the Graph.

        Parameters
        ----------
        styledict : dict, optional
            If provided, allows customization of color and other properties.
        rankdir : str, optional
            The direction of the graph layout. Default is "BT" (bottom to top).
        hide_types : bool, optional
            Whether to hide the types in the visualization. Default is False.
        workflow_view : bool, optional
            Whether to enable the workflow view. Default is False.
        sample_view : bool, optional
            Whether to enable the sample view. Default is False.
        size : tuple, optional
            The size of the visualization. Default is None.
        layout : str, optional
            The name of the layout algorithm for the graph. Default is "neato".

        Returns
        -------
        graphviz.dot.Digraph
            The visualization of the RDF tree.

        Notes
        -----
        The `styledict` parameter allows customization of the visualization style.
        It has the following options:

        BNode:
            color : str
                The color of the BNode boxes.
            shape : str
                The shape of the BNode boxes.
            style : str
                The style of the BNode boxes.

        URIRef:
            color : str
                The color of the URIRef boxes.
            shape : str
                The shape of the URIRef boxes.
            style : str
                The style of the URIRef boxes.

        Literal:
            color : str
                The color of the Literal boxes.
            shape : str
                The shape of the Literal boxes.
            style : str
                The style of the Literal boxes.
        """
        if size is not None:
            size = f"{size[0]},{size[1]}"

        sdict = defstyledict.copy()
        if styledict is not None:
            sdict = _replace_keys(sdict, styledict)

        return visualize_graph(
            self.graph,
            styledict=sdict,
            rankdir=rankdir,
            hide_types=hide_types,
            workflow_view=workflow_view,
            sample_view=sample_view,
            size=size,
            layout=layout,
        )

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

    def close(self, filename, format="json-ld"):
        """
        Close the graph and write to a file

        Parameters
        ----------
        filename: string
            name of output file

        Returns
        -------
        None
        """
        self.write(filename, format=format)

    def archive(self, package_name, format="turtle", compress=True, add_simulations=False):
        """
        Publish a dataset from graph including per atom quantities.

        Parameters:
        -----------
        package_name : str
            The name of the package to be created.
        format : str, optional
            The format in which the dataset should be written. Default is "turtle".
        compress : bool, optional
            Whether to compress the package into a tarball. Default is True.

        Raises:
        -------
        ValueError
            If the package_name already exists or if the tarball already exists.

        Notes:
        ------
        This method creates a package containing a dataset from the graph, including per atom quantities.
        The package consists of a folder named package_name, which contains the dataset and related files.
        If compress is True, the package is compressed into a tarball.

        The method performs the following steps:
        1. Checks if the package_name already exists. If it does, raises a ValueError.
        2. If compress is True, checks if the tarball already exists. If it does, raises a ValueError.
        3. Creates a folder named package_name.
        4. Creates a subfolder named rdf_structure_store within the package folder.
        5. Copies the files associated with each sample to the rdf_structure_store folder, while fixing the paths.
        6. Updates the paths in the graph to point to the copied files.
        7. Writes the dataset to a file named "triples" within the package folder.
        8. If compress is True, compresses the package folder into a tarball.
        9. Removes the package folder.

        """
        # first step make a folder
        if os.path.exists(package_name):
            raise ValueError(f"{package_name} already exists")
        if compress:
            if os.path.exists(f"{package_name}.tar.gz"):
                raise ValueError(f"{package_name} tarball already exists")

        os.mkdir(package_name)
        structure_store = f"{package_name}/rdf_structure_store"
        os.mkdir(structure_store)

        # now go through each sample, and copy the file, at the same time fix the paths
        for sample in self.sample_ids:
            filepath = self.value(URIRef(f"{sample}_Position"), CMSO.hasPath).toPython()
            shutil.copy(filepath, structure_store)

            # now we have to remove the old path, and fix new
            for val in ["Position", "Species"]:
                self.remove((URIRef(f"{sample}_{val}"), CMSO.hasPath, None))

                # assign corrected path
                new_relpath = "/".join(["rdf_structure_store", filepath.split("/")[-1]])
                self.add(
                    (
                        URIRef(f"{sample}_{val}"),
                        CMSO.hasPath,
                        Literal(new_relpath, datatype=XSD.string),
                    )
                )
        #copy simulation files if needed
        if add_simulations:
            sim_store = f"{package_name}/simulation_store"
            os.mkdir(sim_store)
            activities = self.activity_ids
            for activity in activities:
                path = self.value(activity, CMSO.hasPath)
                if path is not None:
                    newpath = "/".join([sim_store, activity.toPython()])
                    shutil.copytree(path, newpath)

                    #remove old path
                    self.remove((activity, CMSO.hasPath, None))
                    self.add((activity, CMSO.hasPath, Literal(newpath, datatype=XSD.string)))

        triple_file = os.path.join(package_name, "triples")
        self.write(triple_file, format=format)

        if compress:
            with tarfile.open(f"{package_name}.tar.gz", "w:gz") as tar:
                tar.add(package_name, arcname=os.path.basename(package_name))
            shutil.rmtree(package_name)

    @classmethod
    def unarchive(
        cls,
        package_name,
        compress=True,
        store="Memory",
        store_file=None,
        identifier="http://default_graph",
        ontology=None,
    ):
        """
        Unarchives a package and returns an instance of the Graph class.

        Parameters
        ----------
        package_name : str
            The name of the package to unarchive.
        compress : bool, optional
            Whether to compress the package. Defaults to True.
        store : str, optional
            The type of store to use. Defaults to "Memory".
        store_file : str, optional
            The file to use for the store. Defaults to None.
        identifier : str, optional
            The identifier for the graph. Defaults to "http://default_graph".
        ontology : str, optional
            The ontology to use. Defaults to None.

        Returns
        -------
        Graph
            An instance of the Graph class.

        Raises
        ------
        FileNotFoundError
            If the package file is not found.
        tarfile.TarError
            If there is an error while extracting the package.

        """
        if compress:
            package_base_name = ".".join(package_name.split(".")[:-2])
            with tarfile.open(package_name) as fin:
                fin.extractall(".")
            # os.remove(package_name)
            # copy things out
        else:
            package_base_name = package_name

        return cls(
            store=store,
            store_file=store_file,
            identifier=identifier,
            graph_file=f"{package_base_name}/triples",
            structure_store=f"{package_base_name}/rdf_structure_store",
            ontology=ontology,
        )

    def query(self, inquery, return_df=True):
        """
        Query the graph using SPARQL

        Parameters
        ----------
        inquery: string
            SPARQL query to be executed

        return_df: bool, optional
            if True, returns the results as a pandas DataFrame. Default is True.
        Returns
        -------
        res: pandas DataFrame
            pandas dataframe results
        """
        res = self.graph.query(inquery)
        if res is not None:
            if return_df:
                for line in inquery.split("\n"):
                    if "SELECT DISTINCT" in line:
                        break
                labels = [x[1:] for x in line.split()[2:]]
                return pd.DataFrame(res, columns=labels)
            else:
                return res
        raise ValueError("SPARQL query returned None")

    def auto_query(
        self,
        source,
        destination,
        return_query=False,
        enforce_types=True,
        return_df=True,
    ):
        """
        Automatically generates and executes a query based on the provided parameters.

        Parameters
        ----------
        source : OntoTerm
            The source of the query.
        destination : OntoTerm
            The destination of the query.
        return_query : bool, optional
            If True, returns the generated query instead of executing it. Defaults to False.
        enforce_types : bool, optional
            If provided, enforces the specified type for the query. Defaults to None.
        return_df: bool, optional
            if True, returns the results as a pandas DataFrame. Default is True.

        Returns
        -------
        pandas DataFrame or str
            The result of the query execution. If `return_query` is True, returns the generated query as a string.
            Otherwise, returns the result of the query execution as a pandas DataFrame.
        """
        query = self.ontology.create_query(
            source, destination, enforce_types=enforce_types
        )
        if return_query:
            return query
        res = self.query(query, return_df=return_df)
        return res

    #################################
    # Methods to interact with sample
    #################################
    def query_sample(
        self, destination, return_query=False, enforce_types=None
    ):
        """
        Query the knowledge graph for atomic scale samples.

        Parameters
        ----------
        destination : OntoTerm
            The destination of the query.
        return_query : bool, optional
            If True, returns the generated query instead of executing it. Defaults to False.
        enforce_types : bool, optional
            If provided, enforces the specified type for the query. Defaults to None.

        Returns
        -------
        pandas DataFrame or str
            The result of the query execution. If `return_query` is True, returns the generated query as a string. Otherwise, returns the result of the query execution as a pandas DataFrame.

        """
        return self.auto_query(
            self.ontology.terms.cmso.AtomicScaleSample,
            destination,
            return_query=return_query,
            enforce_types=enforce_types,
        )

    @property
    def n_samples(self):
        """
        Number of samples in the Graph
        """

        return len([x for x in self.triples((None, RDF.type, CMSO.AtomicScaleSample))])

    @property
    def sample_ids(self):
        """
        Returns a list of all Samples in the graph
        """

        return [x[0] for x in self.triples((None, RDF.type, CMSO.AtomicScaleSample))]

    @property
    def sample_names(self):
        """
        Returns a list of all Sample names in the graph
        """
        samples = [x[0] for x in self.triples((None, RDF.type, CMSO.AtomicScaleSample))]
        samples_names = []
        for sample in samples:
            sample_name = self.value(sample, RDFS.label)
            if sample_name is not None:
                samples_names.append(sample_name.toPython())
            else:
                samples_names.append(sample.toPython())
        return samples_names

    @property
    def samples(self):
        sample_ids = self.sample_ids
        sample_names = self.sample_names
        sample_objects = []
        for ids, name in zip(sample_ids, sample_names):
            sample_objects.append(Sample(name, ids, self))
        return sample_objects

    @property
    def activity_ids(self):
        """
        Returns a list of all Samples in the graph
        """

        return [x[0] for x in self.triples((None, RDF.type, PROV.Activity))]
    
    def iterate_graph(self, item, create_new_graph=False):
        """
        Iterate through the graph starting from the given item.

        Parameters
        ----------
        item : object
            The item to start the iteration from.
        create_new_graph : bool, optional
            If True, create a new KnowledgeGraph object to store the iteration results.
            Default is False. The results are stored in `self.sgraph`.

        Returns
        -------
        None
        """
        if isinstance(item, str):
            item = URIRef(item)

        if create_new_graph:
            self.sgraph = KnowledgeGraph()
        triples = list(self.triples((item, None, None)))
        for triple in triples:
            self.sgraph.graph.add(triple)
            self.iterate_graph(triple[2])

    def get_sample(self, sample, no_atoms=False):
        """
        Get the Sample as a KnowledgeGraph

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
        
        na: int, only retured if no_atoms is True

        """
        if isinstance(sample, str):
            sample = URIRef(sample)

        self.iterate_graph(sample, create_new_graph=True)
        if no_atoms:
            na = self.sgraph.value(sample, CMSO.hasNumberOfAtoms).toPython()
            return self.sgraph, na
        return self.sgraph
    
    def get_label(self, item):
        label = self.graph.value(item, RDFS.label)
        if label is not None:
            return label.toPython()

    def get_sample_label(self, sample):
        label = self.get_label(sample)
        return label
    
    def change_label(self, sample, label):
        self.graph.remove((sample, RDFS.label, None))
        self.graph.add((sample, RDFS.label, Literal(label, datatype=XSD.string)))
    
    def get_system_from_sample(self, sample):
        """
        Get a pyscal :py:class:`atomrdf.structure.System` from the selected sample

        Parameters
        ----------
        sample: string
            sample id

        Returns
        -------
        system: :py:class:`atomrdf.structure.System`
            corresponding system
        """

        simcell = self.value(sample, CMSO.hasSimulationCell)
        cell_vectors = [[], [], []]

        for s in self.triples((simcell, CMSO.hasVector, None)):
            cell_vectors[0].append(self.value(s[2], CMSO.hasComponent_x).toPython())
            cell_vectors[1].append(self.value(s[2], CMSO.hasComponent_y).toPython())
            cell_vectors[2].append(self.value(s[2], CMSO.hasComponent_z).toPython())

        # cell_vectors
        filepath = self.value(URIRef(f"{sample}_Position"), CMSO.hasPath).toPython()
        position_identifier = self.value(
            URIRef(f"{sample}_Position"), CMSO.hasIdentifier
        ).toPython()
        species_identifier = self.value(
            URIRef(f"{sample}_Species"), CMSO.hasIdentifier
        ).toPython()

        # open the file for reading
        with open(filepath, "r") as fin:
            data = json.load(fin)
            positions = data[position_identifier]["value"]
            species = data[species_identifier]["value"]

        atoms = {"positions": positions, "species": species}
        at = Atoms()
        at.from_dict(atoms)
        sys = System()
        sys.box = cell_vectors
        sys.atoms = at
        sys.sample = sample
        sys.graph = self
        sys._name = sample.toPython().split('sample:')[-1]
        return sys

    def to_file(self, sample, filename=None, format="poscar", add_sample_id=True,
        input_data=None, pseudopotentials=None, 
                kspacing=None, kpts=None, 
                koffset=(0, 0, 0), 
                crystal_coordinates=False):
        """
        Save a given sample to a file

        Parameters
        ----------
        sample
            ID of the sample

        filename: string
            name of output file

        format: string, {"lammps-dump","lammps-data", "poscar", 'cif', 'quantum-espresso'}
            or any format supported by ase

        input_data : str, optional
            Additional input data to include in the output file. Defaults to None.
            Only valid for quantum-espresso format. See ASE write docs for more information.
        pseudopotentials : str, optional
            The path to the pseudopotentials file. Defaults to None.
            Only valid for quantum-espresso format. See ASE write docs for more information.
        kspacing : float, optional
            The k-spacing value to include in the output file. Defaults to None.
            Only valid for quantum-espresso format. See ASE write docs for more information.
        kpts : list, optional
            A list of k-points to include in the output file. Defaults to None.
            Only valid for quantum-espresso format. See ASE write docs for more information.
        koffset : tuple, optional
            The k-offset values to include in the output file. Defaults to (0, 0, 0).
            Only valid for quantum-espresso format. See ASE write docs for more information.
        crystal_coordinates : bool, optional
            Whether to include crystal coordinates in the output file. Defaults to False.
            Only valid for quantum-espresso format. See ASE write docs for more information.

        Returns
        -------
        None
        """

        if filename is None:
            filename = os.path.join(os.getcwd(), "out")

        sys = self.get_system_from_sample(sample)
        sys.to_file(filename=filename, format=format, add_sample_id=add_sample_id, input_data=input_data, 
                pseudopotentials=pseudopotentials, kspacing=kspacing, 
                kpts=kpts, koffset=koffset, crystal_coordinates=crystal_coordinates)

    def enable_workflow(self, workflow_object, workflow_environment=None, workflow_module=None):
        self.workflow.inform_graph(workflow_object, 
                        workflow_environment=workflow_environment, 
                        workflow_module=workflow_module)
        
    def add_workflow(self, job, workflow_environment=None, workflow_module=None, job_dicts=None,
                    add_intermediate_jobs=False):
        self.workflow.to_graph(job, workflow_environment=workflow_environment, 
                            workflow_module=workflow_module, 
                            job_dicts=job_dicts,
                            add_intermediate_jobs=add_intermediate_jobs)
    
    def find_property(self, label):
        prop_list = list(self.graph.triples((None, RDFS.label, label)))
        if len(prop_list) == 0:
            raise RuntimeError(f'Property {label} not found in the graph')
        prop = prop_list[0][0]
        return prop        

    def get_string_label(self, item):
        label = self.get_label(item)

        if label is None:
            try:
                label = str(item.toPython())
            except:
                label = str(item)
        
        if "activity" in label:
            method = self.value(item, ASMO.hasComputationalMethod)
            if method is not None:
                method_name = self.value(method, RDF.type)
                if method_name is not None:
                    label = method_name.toPython().split("/")[-1]
        return label

    def _add_to_dict(self, prop, indict):
        name = _name(prop)
        if name not in indict.keys():
            indict[name] = {}
            indict[name]['found'] = False
            indict[name]['label'] = self.get_string_label(prop)

    def _get_ancestor(self, prop, prov):
        #note that only one operation and parent are present!
        if isinstance(prop, str):
            prop = URIRef(prop)
        propname = _name(prop)

        operation = [x[1] for x in self.triples((prop, ASMO.wasCalculatedBy, None))]
        if len(operation) > 0:
            parent = [x[2] for x in self.triples((prop, ASMO.wasCalculatedBy, None))]
            operation = operation[0]
            parent = parent[0]        
            prov[propname]['operation'] = 'output_parameter'
            prov[propname]['inputs'] = {}
            prov[propname]['inputs']['0'] = _name(parent)
            self._add_to_dict(parent, prov)
            prov[_name(parent)]['inputs'] = {}
            associated_samples = [x[0] for x in self.triples((None, PROV.wasGeneratedBy, parent))]
            for count, sample in enumerate(associated_samples):
                prov[_name(parent)]['inputs'][str(count)] = _name(sample)
                self._add_to_dict(sample, prov)
            prov[_name(parent)]['found'] = True
            prov[_name(parent)]['operation'] = 'sample_for_activity'
                
        else:
            operation = [x[1] for x in self.triples((None, None, prop))]
            parent = [list(self.triples((None, op, prop)))[0][0] for op in operation]
            if len(operation) == 0:
                prov[propname]['found'] = True
                return prov
            operation = operation[0]
            parent = parent[0]

            if operation.toPython() == "http://purls.helmholtz-metadaten.de/asmo/hasInputParameter":
                #print(f'we ran 1 for {operation.toPython()}')
                prov[propname]['operation'] = 'input_parameter'
                prov[propname]['inputs'] = {}
                prov[propname]['inputs']['0'] = _name(parent)
                self._add_to_dict(parent, prov)
                prov[_name(parent)]['inputs'] = {}
                associated_samples = [x[0] for x in self.triples((None, PROV.wasGeneratedBy, parent))]
                for count, sample in enumerate(associated_samples):
                    prov[_name(parent)]['inputs'][str(count)] = _name(sample)
                    self._add_to_dict(sample, prov)
                prov[_name(parent)]['found'] = True
                prov[_name(parent)]['operation'] = 'sample_for_activity'

            elif operation.toPython() == "http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty":
                #print(f'we ran 2 for {operation.toPython()}')
                prov[propname]['operation'] = 'output_parameter'
                prov[propname]['inputs'] = {}
                prov[propname]['inputs']['0'] = _name(parent)
                self._add_to_dict(parent, prov)            
                prov[_name(parent)]['found'] = True
                prov[_name(parent)]['operation'] = 'sample_output'
            
            elif operation == MATH.hasSum:
                addends = list(x[2] for x in self.triples((parent, MATH.hasAddend, None)))
                prov[propname]['operation'] = 'addition'
                prov[propname]['inputs'] = {}
                for count, term in enumerate(addends):
                    prov[propname]['inputs'][f'{count}'] = _name(term)
                    self._add_to_dict(term, prov)
            
            elif operation == MATH.hasDifference:
                minuend = self.value(parent, MATH.hasMinuend)
                subtrahend = self.value(parent, MATH.hasSubtrahend)
                prov[propname]['operation'] = 'subtraction'
                prov[propname]['inputs'] = {}
                prov[propname]['inputs']['0'] = _name(minuend)
                prov[propname]['inputs']['1'] = _name(subtrahend)
                self._add_to_dict(minuend, prov)
                self._add_to_dict(subtrahend, prov)
            
            elif operation == MATH.hasProduct:
                factors = list(x[2] for x in self.triples((parent, MATH.hasFactor, None)))
                prov[propname]['operation'] = 'multiplication'
                prov[propname]['inputs'] = {}
                for count, term in enumerate(factors):
                    prov[propname]['inputs'][f'{count}'] = _name(term)
                    self._add_to_dict(term, prov)
            
            elif operation == MATH.hasQuotient:
                divisor = self.value(parent, MATH.hasDivisor)
                dividend = self.value(parent, MATH.hasDividend)
                prov[propname]['operation'] = 'division'
                prov[propname]['inputs'] = {}
                prov[propname]['inputs']['0'] = _name(divisor)
                prov[propname]['inputs']['1'] = _name(dividend)
                self._add_to_dict(divisor, prov)
                self._add_to_dict(dividend, prov)

        prov[propname]['found'] = True
        return prov
    
    def generate_provenance(self, prop=None, label=None, visualize=False):
        if (prop is None) and (label is None):
            raise ValueError('Either prop or label must be provided')
        
        if prop is None:
            prop = self.find_property(label)
        
        name = _name(prop)
        prov = {}
        self._add_to_dict(prop, prov)

        done = False
        while not done:
            done = True
            keys = list(prov.keys()) 
            for prop in keys:
                if not prov[prop]['found']:
                    prov = self._get_ancestor(prop, prov)
                    done = False
        
        if visualize:
            return visualize_provenance(prov)
        
        return prov
    

    def query_structure_from_mp(self, api_key, chemical_system=None, material_ids=None, is_stable=True,
                        conventional=True,
                        add_to_graph=True):
        docs = amp.query_mp(api_key, 
                            chemical_system = chemical_system, 
                            material_ids = material_ids, 
                            is_stable = is_stable)
        structures = []
        for doc in docs:
            struct = doc['structure']
            if conventional:
                aseatoms = struct.to_conventional().to_ase_atoms()
            else:
                aseatoms = struct.to_primitive().to_ase_atoms()
            
            sys = System.read.ase(aseatoms)
            
            symmetry = doc['symmetry']

            if add_to_graph:
                sys.graph = self
                sys.to_graph()
                
                targets = [None, symmetry['symbol'], symmetry['number'], None, None, None]
                sys._add_crystal_structure(targets=targets)

                #add energy
                self.add_calculated_quantity(sys.sample, 
                                            'EnergyPerAtom', 
                                            doc['energy_per_atom'], 
                                            unit='EV')  
                structures.append(sys)
        if len(structures) == 1:
            return structures[0]
        else:
            return structures                     
