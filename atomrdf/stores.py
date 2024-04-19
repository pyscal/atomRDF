from rdflib.store import NO_STORE, VALID_STORE
from rdflib import plugin
from rdflib import Graph, Literal

import os

# special methods; for supporting workflow envs
from atomrdf.workflow import inform_graph


def create_store(kg, store, identifier, store_file=None, structure_store=None):
    """
    Create a store based on the given parameters.

    Parameters:
    -----------
    kg : KnowledgeGraph
        The knowledge graph object.
    store : str or Project
        The type of store to create. It can be either "Memory", "SQLAlchemy", or a pyiron Project object.
    identifier : str
        The identifier for the store.
    store_file : str, optional
        The file path to store the data (only applicable for certain store types).
    structure_store : str, optional
        The structure store to use (only applicable for certain store types).

    Raises:
    -------
    ValueError
        If an unknown store type is provided.

    """
    kg.store_file = store_file
    if store == "Memory":
        store_memory(
            kg,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )
    elif store == "SQLAlchemy":
        store_alchemy(
            kg,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )
    elif type(store).__name__ == "Project":
        store_pyiron(
            kg,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )
    else:
        raise ValueError("Unknown store found!")


def store_memory(kg, store, identifier, store_file=None, structure_store=None):
    """
    Store the knowledge graph in memory.

    Parameters
    ----------
    kg : KnowledgeGraph
        The knowledge graph to be stored.
    store : str
        The type of store to use for storing the graph.
    identifier : str
        The identifier for the graph.
    store_file : str, optional
        The file to store the graph in. Defaults to None.
    structure_store : str, optional
        The structure store to use. Defaults to None.

    Returns
    -------
    None
    """
    graph = Graph(store="Memory", identifier=identifier)
    kg.graph = graph
    kg.structure_store = _setup_structure_store(structure_store=structure_store)


def store_alchemy(kg, store, identifier, store_file=None, structure_store=None):
    """
    Store the knowledge graph using SQLAlchemy.

    Parameters
    ----------
    kg : KnowledgeGraph
        The knowledge graph to be stored.
    store : str
        The type of store to be used.
    identifier : str
        The identifier for the graph.
    store_file : str, optional
        The file path for the store. Required if store is not 'memory'.
    structure_store : str, optional
        The structure store to be used.

    Raises
    ------
    ValueError
        If store_file is None and store is not 'memory'.

    Returns
    -------
    None
    """
    _check_if_sqlalchemy_is_available()
    if store_file is None:
        raise ValueError("store file is needed if store is not memory")

    kg.graph = Graph(store="SQLAlchemy", identifier=identifier)
    uri = Literal(f"sqlite:///{store_file}")
    kg.graph.open(uri, create=True)
    kg.structure_store = _setup_structure_store(structure_store=structure_store)


def store_pyiron(kg, store, identifier, store_file=None, structure_store=None):
    """
    Store the pyiron knowledge graph in a database.

    Parameters
    ----------
    kg : pyiron.atomistics.structure.pyiron_atomistics.structure.AtomisticStructure
        The pyiron knowledge graph to be stored.
    store : pyiron.atomistics.structure.pyiron_atomistics.structure.AtomisticStructure
        The store object where the knowledge graph will be stored.
    identifier : str
        The identifier for the knowledge graph.
    store_file : str, optional
        The path to the store file. If not provided, a default path will be used.
    structure_store : str, optional
        The path to the structure store. If not provided, a default path will be used.

    Returns
    -------
    None

    """
    structure_store = os.path.join(store.path, "rdf_structure_store")
    kg.structure_store = _setup_structure_store(structure_store=structure_store)
    store_file = os.path.join(store.path, f"{store.name}.db")
    store_alchemy(kg, store, identifier, store_file, structure_store=structure_store)
    # finally update project object
    inform_graph(store, kg)


def _check_if_sqlalchemy_is_available():
    try:
        import sqlalchemy as sa
    except ImportError:
        raise RuntimeError("Please install the sqlalchemy package")
    try:
        import rdflib_sqlalchemy as rsa
    except ImportError:
        raise RuntimeError(
            "Please install the rdllib-sqlalchemy package. The development version is needed, please do pip install git+https://github.com/RDFLib/rdflib-sqlalchemy.git@develop"
        )


def _setup_structure_store(structure_store=None):
    if structure_store is None:
        structure_store = os.path.join(os.getcwd(), "rdf_structure_store")
    if not os.path.exists(structure_store):
        os.mkdir(structure_store)
    return structure_store
