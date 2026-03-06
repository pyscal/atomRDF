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

from pyscal3.atoms import Atoms

from atomrdf.visualize import visualize_graph, visualize_provenance
from atomrdf.ontology import read_ontology

import atomrdf.properties as prp
from atomrdf.stores import create_store, purge
import atomrdf.json_io as json_io
import atomrdf.mp as amp


from atomrdf.namespace import (
    Namespace,
    CMSO,
    PLDO,
    PODO,
    ASMO,
    PROV,
    MATH,
    CDCO,
    UNSAFECMSO,
    UNSAFEASMO,
    Literal,
)

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
    input_string = re.sub(r"\W", "_", input_string)
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
    add(triple, validate=True)
        Add a triple to the knowledge graph.
    triples(triple)
        Return the triples in the knowledge graph that match the given triple pattern.
    query(source, destinations=None, return_df=True, num_paths=1, limit=None)
        Execute a SPARQL query on the knowledge graph using tools4RDF.
    get_sample_as_structure(sample_id)
        Retrieve a sample from the graph as an AtomicScaleSample object.
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

        self._store = store
        self._identifier = identifier
        self._store_file = store_file

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

    def purge(self, force=False):
        """
        Remove all information from the KnowledgeGraph.

        Parameters
        ----------
        force : bool, optional
            Whether to proceed with purging the graph. Default is False.

        Returns
        -------
        None

        Notes
        -----
        This method removes all information from the KnowledgeGraph. If the `force` parameter is set to False, a warning is issued before proceeding with the purging.
        """
        if not force:
            warnings.warn(
                "This will remove all information from the KnowledgeGraph. Call with force=True to proceed."
            )
            return
        else:
            # Clean up structure store files referenced by this graph
            # Query for all files referenced via CMSO.hasPath
            file_paths = set()
            for triple in self.triples((None, CMSO.hasPath, None)):
                filepath = triple[2].toPython()
                if filepath:
                    full_path = os.path.join(
                        self.structure_store, os.path.basename(filepath)
                    )
                    if os.path.exists(full_path):
                        file_paths.add(full_path)

            # Remove the files
            for filepath in file_paths:
                os.remove(filepath)

            # Close the current store before destroying and recreating it
            # (required for file-backed stores like Oxigraph that use a lock file)
            try:
                self.graph.close()
            except Exception:
                pass

            graph = purge(self._store, self._identifier, self._store_file)
            self.graph = graph
            self._n_triples = 0

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

    def _is_uriref(self, term):
        return type(term).__name__ == "URIRef"

    def _is_bnode(self, term):
        return not term.toPython().startswith("http")

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

        # if destination_range == "string":
        #    destination_range = "str"
        # elif destination_range == "integer":
        #    destination_range = "int"

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

    def objects(self, arg1, arg2):
        modified_double = self._modify_triple((arg1, arg2))
        return self.graph.objects(modified_double[0], modified_double[1])

    def query(self, source, destinations=None, return_df=True, num_paths=1, limit=None):
        """
        Execute a SPARQL query on the knowledge graph.

        This method supports two query modes:
        1. Raw SPARQL query strings (passed as source parameter)
        2. Ontology-based queries using tools4RDF (source as OntoTerm)

        Parameters
        ----------
        source : str or OntoTerm
            If str: Raw SPARQL query string to execute directly.
            If OntoTerm: The source ontology term from which paths are to be queried.
            Access terms via self.ontology.terms (e.g., self.ontology.terms.cmso.AtomicScaleSample).
        destinations : list of OntoTerm or OntoTerm, optional
            One or more destination ontology terms to which paths are to be queried.
            Can be a single term or a list of terms. If None, all properties of the source are returned.
            Only used when source is an OntoTerm.
        return_df : bool, default=True
            If True, returns results as a pandas DataFrame. Otherwise, returns raw query results.
        num_paths : int, default=1
            The number of paths to retrieve for each query when multiple paths exist.
            Only used when source is an OntoTerm.
        limit : int, optional
            The maximum number of results to return. If None, no limit is applied.
            Only used when source is an OntoTerm.

        Returns
        -------
        pandas.DataFrame or list or None
            If return_df is True, returns a pandas DataFrame with query results.
            If return_df is False, returns a list of query results.
            Returns None if no results are found.

        Examples
        --------
        Query with raw SPARQL string:

        >>> query = '''
        ... PREFIX cmso: <http://purls.helmholtz-metadaten.de/cmso/>
        ... SELECT DISTINCT ?symbol
        ... WHERE {
        ...     ?sample cmso:hasNumberOfAtoms ?number .
        ...     ?sample cmso:hasMaterial ?material .
        ...     ?material cmso:hasStructure ?structure .
        ...     ?structure cmso:hasSpaceGroupSymbol ?symbol .
        ... FILTER (?number="4"^^xsd:integer)
        ... }'''
        >>> df = kg.query(query)

        Query for all AtomicScaleSamples with their space group symbols:

        >>> kg = KnowledgeGraph()
        >>> df = kg.query(
        ...     kg.ontology.terms.cmso.AtomicScaleSample,
        ...     [kg.ontology.terms.cmso.hasSpaceGroupSymbol]
        ... )

        Query with filters (using == operator on terms):

        >>> df = kg.query(
        ...     kg.ontology.terms.cmso.AtomicScaleSample,
        ...     [kg.ontology.terms.cmso.hasNumberOfAtoms == 4]
        ... )

        Notes
        -----
        When using ontology terms, this method uses tools4RDF to automatically generate
        SPARQL queries based on the ontology structure. It handles namespace management,
        path finding between ontology terms, and result formatting automatically.
        """
        # Check if source is a raw SPARQL query string
        if isinstance(source, str):
            res = self.graph.query(source)
            if res is not None and return_df:
                # Extract column names from SELECT clause
                import re

                select_match = re.search(
                    r"SELECT\s+(?:DISTINCT\s+)?(.+?)\s+WHERE",
                    source,
                    re.IGNORECASE | re.DOTALL,
                )
                if select_match:
                    # Extract variable names (anything starting with ?)
                    variables = re.findall(r"\?(\w+)", select_match.group(1))
                    return pd.DataFrame(res, columns=variables)
                else:
                    # Fallback: return as DataFrame without column names
                    return pd.DataFrame(res)
            return res
        else:
            # Use tools4RDF for ontology-based queries
            return self.ontology.query(
                self.graph,
                source,
                destinations=destinations,
                return_df=return_df,
                num_paths=num_paths,
                limit=limit,
            )

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

    def close_store(self):
        """
        Release the underlying store (close file handles and locks).

        This is a no-op for the in-memory store. For file-backed stores
        (Oxigraph, SQLAlchemy) it releases the file lock so the same store
        directory can be reopened in the same process or by another process.

        Returns
        -------
        None
        """
        self.graph.close()

    def __del__(self):
        try:
            self.graph.close()
        except Exception:
            pass

    def archive(
        self, package_name, format="turtle", compress=True, add_simulations=False
    ):
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
        copied_basenames = set()

        # First, attempt to copy the canonical sample files as before (best-effort)
        for sample in self.sample_ids:
            filepath = self.value(URIRef(f"{sample}_Position"), CMSO.hasPath)
            if filepath is None:
                continue
            filepath = filepath.toPython()
            # filepath has to be fixed with the correct prefix as needed
            srcpath = os.path.join(self.structure_store, os.path.basename(filepath))
            if os.path.exists(srcpath):
                if os.path.basename(filepath) not in copied_basenames:
                    shutil.copy(srcpath, structure_store)
                    copied_basenames.add(os.path.basename(filepath))

            # now we have to remove the old path, and fix new for Position/Species nodes
            for val in ["Position", "Species"]:
                self.remove((URIRef(f"{sample}_{val}"), CMSO.hasPath, None))

                # assign corrected path
                new_relpath = "/".join(
                    ["rdf_structure_store", os.path.basename(filepath)]
                )
                self.add(
                    (
                        URIRef(f"{sample}_{val}"),
                        CMSO.hasPath,
                        Literal(new_relpath, datatype=XSD.string),
                    )
                )

        # Additionally, copy any other files referenced by CMSO.hasPath (e.g. calculated-property files)
        for subj, pred, obj in list(self.triples((None, CMSO.hasPath, None))):
            try:
                path = obj.toPython()
            except Exception:
                continue
            basename = os.path.basename(path)
            src = os.path.join(self.structure_store, basename)
            if not os.path.exists(src):
                # nothing to copy for this path
                continue
            if basename in copied_basenames:
                # already copied
                # still ensure triple points to rdf_structure_store
                self.remove((subj, CMSO.hasPath, None))
                new_relpath = "/".join(["rdf_structure_store", basename])
                self.add(
                    (subj, CMSO.hasPath, Literal(new_relpath, datatype=XSD.string))
                )
                continue

            shutil.copy(src, structure_store)
            copied_basenames.add(basename)

            # replace path triple with corrected relative path inside package
            self.remove((subj, CMSO.hasPath, None))
            new_relpath = "/".join(["rdf_structure_store", basename])
            self.add((subj, CMSO.hasPath, Literal(new_relpath, datatype=XSD.string)))
        # copy simulation files if needed
        if add_simulations:
            sim_store = f"{package_name}/simulation_store"
            os.mkdir(sim_store)
            activities = self.activity_ids
            for activity in activities:
                path = self.value(activity, CMSO.hasPath)
                if path is not None:
                    newpath = "/".join([sim_store, activity.toPython()])
                    shutil.copytree(path, newpath)

                    # remove old path
                    self.remove((activity, CMSO.hasPath, None))
                    self.add(
                        (activity, CMSO.hasPath, Literal(newpath, datatype=XSD.string))
                    )

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
            with tarfile.open(package_name) as fin:
                # Derive the top-level directory from the archive members
                # instead of string-splitting the path (which breaks with
                # absolute paths like /full/path/dataset.tar.gz).
                top_dirs = {m.name.split("/")[0] for m in fin.getmembers()}
                if len(top_dirs) == 1:
                    package_base_name = top_dirs.pop()
                else:
                    # Fallback: strip extensions from the basename
                    package_base_name = os.path.basename(package_name)
                    for ext in (".tar.gz", ".tar.bz2", ".tar"):
                        if package_base_name.endswith(ext):
                            package_base_name = package_base_name[: -len(ext)]
                            break
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

    def merge_archive(self, package_name, compress=True, format="turtle"):
        """
        Merge an archived dataset into this KnowledgeGraph.

        Unlike ``unarchive`` (which creates a new graph), this method loads
        the triples and structure-store files from an existing archive into
        the *current* graph so that multiple datasets can be combined
        incrementally::

            kg = KnowledgeGraph()
            kg.merge_archive("dataset_1_GB.tar.gz")
            kg.merge_archive("dataset_2_GB.tar.gz")
            # kg now contains both datasets

        Parameters
        ----------
        package_name : str
            Path to the archive.  When *compress* is True (default) this
            should be a ``.tar.gz`` file; otherwise the name of an already-
            extracted directory.
        compress : bool, optional
            Whether *package_name* is a compressed tarball.  Default True.
        format : str, optional
            RDF serialisation format of the ``triples`` file inside the
            archive.  Default ``"turtle"``.

        Notes
        -----
        * Structure-store JSON files from the archive are copied into
          ``self.structure_store``.  UUID-based filenames make collisions
          extremely unlikely; a warning is emitted if a file already exists
          and it is silently skipped (the existing copy wins).
        * After parsing, every ``CMSO.hasPath`` triple that still references
          the archive-internal ``rdf_structure_store/`` prefix is rewritten
          to point at ``self.structure_store``.
        """
        # --- 1. Extract if compressed ----------------------------------------
        if compress:
            with tarfile.open(package_name) as fin:
                # The archive's top-level directory name is used as package_base.
                # extractall(".") puts it relative to cwd, so we must derive
                # package_base from the archive member names rather than the
                # (possibly absolute) package_name path.
                top_dirs = {m.name.split("/")[0] for m in fin.getmembers()}
                fin.extractall(".")
            if len(top_dirs) == 1:
                package_base = top_dirs.pop()
            else:
                # fallback: strip .tar.gz from basename
                package_base = os.path.basename(package_name)
                for suffix in (".tar.gz", ".tgz"):
                    if package_base.endswith(suffix):
                        package_base = package_base[: -len(suffix)]
                        break
        else:
            package_base = package_name

        triple_file = os.path.join(package_base, "triples")
        archive_ss = os.path.join(package_base, "rdf_structure_store")

        # --- 2. Copy structure-store files -----------------------------------
        if os.path.isdir(archive_ss):
            for fname in os.listdir(archive_ss):
                src = os.path.join(archive_ss, fname)
                dst = os.path.join(self.structure_store, fname)
                if os.path.exists(dst):
                    warnings.warn(
                        f"merge_archive: '{fname}' already exists in "
                        f"structure store — skipping (existing copy kept)",
                        UserWarning,
                        stacklevel=2,
                    )
                    continue
                shutil.copy(src, dst)

        # --- 3. Parse triples ------------------------------------------------
        if os.path.exists(triple_file):
            self.graph.parse(triple_file, format=format)

        # --- 4. Rewrite CMSO.hasPath to point at self.structure_store --------
        archive_prefix = "rdf_structure_store/"
        for subj, pred, obj in list(self.triples((None, CMSO.hasPath, None))):
            try:
                path = obj.toPython()
            except Exception:
                continue
            if not isinstance(path, str):
                continue
            basename = os.path.basename(path)
            # Rewrite any path whose basename exists in self.structure_store
            resolved = os.path.join(self.structure_store, basename)
            if os.path.exists(resolved):
                new_relpath = os.path.relpath(resolved, os.getcwd())
                if path != new_relpath:
                    self.remove((subj, CMSO.hasPath, obj))
                    self.add(
                        (
                            subj,
                            CMSO.hasPath,
                            Literal(new_relpath, datatype=XSD.string),
                        )
                    )

        # --- 5. Clean up extracted folder ------------------------------------
        if compress and os.path.isdir(package_base):
            shutil.rmtree(package_base)

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
        Returns a list of all Sample names in the graph.
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

    def _is_of_type(self, item, target_item):
        """
        Check if an item is of a specific type

        item - direct from graph, comparison only makes sense if it is URIRef
            if item is node with https - direct comparison
            if not - check the type of the item

        target_item - URIRef or OntoTerm
        """
        if not self._is_uriref(item):
            return False

        if self._is_bnode(item):
            rdftype = self.value(item, RDF.type)
            if rdftype is not None:
                rdftype = rdftype.toPython()
        else:
            rdftype = item.toPython()

        target_type = target_item.toPython()
        return rdftype == target_type

    def get_sample_as_structure(self, sample_id):
        """
        Retrieve a sample from the graph as an AtomicScaleSample object.

        Parameters
        ----------
        sample_id : str or URIRef
            The ID of the sample to retrieve

        Returns
        -------
        AtomicScaleSample
            The sample as an AtomicScaleSample pydantic object

        Examples
        --------
        >>> kg = KnowledgeGraph()
        >>> sample = kg.get_sample_as_structure('sample:123')
        >>> atoms = sample.to_structure()  # Convert to ASE Atoms
        >>> sample.to_file('output.lmp', format='lammps-dump')
        """
        from atomrdf.datamodels.structure import AtomicScaleSample

        if isinstance(sample_id, str):
            sample_id = (
                sample_id if sample_id.startswith("sample:") else f"sample:{sample_id}"
            )

        return AtomicScaleSample.from_graph(self, sample_id)

    def to_file(
        self,
        sample,
        filename,
        format="lammps-data",
        copy_from=None,
        pseudo_files=None,
    ):
        """
        Write a sample structure to a file.

        Parameters
        ----------
        sample : str or URIRef
            Sample ID

        filename : str
            Name of the output file

        format : str, optional
            Format of the output file. Default is 'lammps-data'.
            Any format supported by ASE can be used.

        copy_from : str, optional
            If provided, input options for quantum-espresso format will be copied from
            the given file. Structure specific information will be replaced.
            Note that the validity of input file is not checked.

        pseudo_files : list, optional
            If provided, add the pseudopotential filenames to file.
            Should be in alphabetical order of chemical species symbols.

        Returns
        -------
        None

        Examples
        --------
        >>> kg = KnowledgeGraph()
        >>> kg.to_file('sample:123', 'output.lmp', 'lammps-data')
        >>> kg.to_file('sample:456', 'POSCAR', 'vasp')
        """
        sample_obj = self.get_sample_as_structure(sample)
        sample_obj.to_file(
            outfile=filename,
            format=format,
            copy_from=copy_from,
            pseudo_files=pseudo_files,
        )

    def get_label(self, item):
        label = self.graph.value(item, RDFS.label)
        if label is not None:
            return label.toPython()

    def get_string_label(self, item):
        label = self.get_label(item)

        if label is None:
            try:
                label = str(item.toPython())
            except:
                label = str(item)

        if "simulation" in label:
            method = self.value(item, ASMO.hasComputationalMethod)
            if method is not None:
                method_name = self.value(method, RDF.type)
                if method_name is not None:
                    label = method_name.toPython().split("/")[-1]
        return label
