from rdflib.store import NO_STORE, VALID_STORE
from rdflib import plugin
from rdflib import Graph
from atomrdf.namespace import Literal

import os
import shutil

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
    if store in ["Memory", "memory"]:
        store_memory(
            kg,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )
    elif store in ["SQLAlchemy", "db", "database", "sqlalchemy"]:
        store_alchemy(
            kg,
            store,
            identifier,
            store_file=store_file,
            structure_store=structure_store,
        )
    elif store in ["Oxigraph", "oxigraph"]:
        store_oxigraph(
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

def store_oxigraph(kg, store, identifier, store_file=None, structure_store=None):
    """
    Store the knowledge graph using Oxigraph (via oxrdflib).

    Parameters
    ----------
    kg : KnowledgeGraph
        The knowledge graph to be stored.
    store : str
        The type of store to be used.
    identifier : str or URIRef
        The URI identifier for the named graph. Must be consistent across
        open/reopen calls to retrieve the same triples.
    store_file : str, optional
        Directory path for the persistent on-disk Oxigraph store.
        If None, an in-memory store is used (data is lost when the object
        is garbage-collected).
    structure_store : str, optional
        The structure store to be used.

    Raises
    ------
    RuntimeError
        If oxrdflib is not installed.

    Returns
    -------
    None
    """
    _check_if_oxrdflib_is_available()
    from rdflib import Graph, URIRef as _URIRef

    if isinstance(identifier, str):
        identifier = _URIRef(identifier)

    graph = Graph(store="Oxigraph", identifier=identifier)

    if store_file is not None:
        create = not os.path.exists(store_file)
        graph.open(store_file, create=create)

    # Oxigraph strictly validates prefix IRIs passed to its SPARQL engine.
    # atomRDF binds some prefixes to local file paths (e.g. cmso.owl) which
    # are not valid absolute IRIs.  Patch the store's bind() to silently
    # skip any namespace whose IRI has no URI scheme.
    _orig_bind = graph.store.bind.__func__  # type: ignore[attr-defined]

    def _safe_bind(self_store, prefix, namespace, override=True):
        ns_str = str(namespace)
        # Skip namespaces that are bare file paths (no scheme)
        if "://" not in ns_str and not ns_str.startswith("urn:"):
            return
        _orig_bind(self_store, prefix, namespace, override)

    import types
    graph.store.bind = types.MethodType(_safe_bind, graph.store)

    kg.graph = graph
    kg.structure_store = _setup_structure_store(structure_store=structure_store)


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


def _check_if_oxrdflib_is_available():
    try:
        import oxrdflib  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "Please install the oxrdflib package: pip install oxrdflib"
        )


def _setup_structure_store(structure_store=None):
    if structure_store is None:
        structure_store = os.path.join(os.getcwd(), "rdf_structure_store")
    if not os.path.exists(structure_store):
        os.mkdir(structure_store)
    return structure_store

def purge(store, identifier, store_file):
    if store in ["Memory", "memory"]:
        return _purge_memory(identifier, store_file)

    elif store in ["SQLAlchemy", "db", "database", "sqlalchemy"]:
        return _purge_alchemy(identifier, store_file)

    elif store in ["Oxigraph", "oxigraph"]:
        return _purge_oxigraph(identifier, store_file)

    else:
        raise ValueError("Unknown store found!")    


def _purge_memory(identifier, store_file):
    graph = Graph(store="Memory", identifier=identifier)
    return graph

def _purge_alchemy(identifier, store_file):
    os.remove(store_file)
    graph = Graph(store="SQLAlchemy", identifier=identifier)
    uri = Literal(f"sqlite:///{store_file}")
    graph.open(uri, create=True)
    return graph


def _purge_oxigraph(identifier, store_file):
    _check_if_oxrdflib_is_available()
    from rdflib import Graph, URIRef as _URIRef

    if isinstance(identifier, str):
        identifier = _URIRef(identifier)

    if store_file is not None and os.path.exists(store_file):
        shutil.rmtree(store_file)

    graph = Graph(store="Oxigraph", identifier=identifier)
    if store_file is not None:
        graph.open(store_file, create=True)
    return graph