import graphviz
import os
from rdflib import BNode, URIRef, Namespace, Literal
import uuid
import json


def get_title_from_BNode(x):
    return x.toPython()


def get_string_from_URI(x):
    """
    Extract a presentable string from URI.

    Parameters
    ----------
    x : rdflib.term.URIRef
        The URI object to extract the string from.

    Returns
    -------
    tuple
        A tuple containing the presentable string representation of the URI and its type.
        The string representation is the last part of the URI after splitting by '#' or '/'.
        The type can be either "URIRef" or "BNode".
    """
    raw = x.toPython()
    # first try splitting by #
    rawsplit = raw.split("#")
    if len(rawsplit) > 1:
        return rawsplit[-1], "URIRef"

    # try splitting by = for chebi values
    if "CHEBI" in raw:
        rawsplit = raw.split("=")
        rawsplit = rawsplit[-1].split(":")
        if len(rawsplit) > 1:
            return ".".join(rawsplit[-2:]), "URIRef"

    if "sample:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "activity:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "operation:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "property:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    # just a normal url split now
    rawsplit = raw.split("/")
    if len(rawsplit) > 1:
        return ".".join(rawsplit[-2:]), "URIRef"

    # none of the conditions worked, which means it's a hex string
    return raw, "BNode"


def parse_object(x):
    """
    Parse the given object and return its title and type.

    Parameters
    ----------
    x : RDF term
        The RDF term to parse.

    Returns
    -------
    tuple
        A tuple containing the title of the object and its type.

    """
    if isinstance(x, BNode):
        return get_title_from_BNode(x), "BNode"
    elif isinstance(x, URIRef):
        return get_string_from_URI(x)
    elif isinstance(x, Literal):
        return str(x.title()), "Literal"


styledict = {
    "BNode": {"color": "#ffe6ff", "shape": "box", "style": "filled"},
    "URIRef": {"color": "#ffffcc", "shape": "box", "style": "filled"},
    "Literal": {"color": "#e6ffcc", "shape": "ellipse", "style": "filled"},
}


def _switch_box(box):
    if box == "box":
        return "rectangle"
    # remember that only boxes will be used, circles no!


def _fix_id(string1, istype1):
    if istype1 == "Literal":
        id1 = str(uuid.uuid4())
    else:
        id1 = string1
    return id1


def visualize_graph(
    g,
    styledict=styledict,
    rankdir="TB",
    hide_types=False,
    workflow_view=False,
    sample_view=False,
    size=None,
    layout="dot",
):
    """
    Visualizes a graph using Graphviz.

    Parameters
    ----------
    g : dict
        The graph to visualize.
    styledict : dict, optional
        A dictionary containing styles for different types of nodes and edges. Default is `styledict`.
    rankdir : str, optional
        The direction of the graph layout. Default is "TB" (top to bottom).
    hide_types : bool, optional
        Whether to hide nodes with the "type" attribute. Default is False.
    workflow_view : bool, optional
        Whether to enable the workflow view. Default is False.
    sample_view : bool, optional
        Whether to enable the sample view. Default is False.
    size : str, optional
        The size of the graph. Default is None.
    layout : str, optional
        The layout algorithm to use. Default is "dot".

    Returns
    -------
    dot : graphviz.Digraph
        The graph visualization.
    """
    dot = graphviz.Digraph()

    dot.attr(
        rankdir=rankdir,
        style="filled",
        size=size,
        layout=layout,
        overlap="false",
    )

    for k in g:
        string1, istype1 = parse_object(k[0])
        string2, istype2 = parse_object(k[2])
        string3, istype = parse_object(k[1])

        plot = True

        if workflow_view:
            # we collapse sample information
            # if cmso.connector is found, only use it is it is cmso.hasCalculated
            # all sub sample props, indicated by sample_x_jsjsj will be ignored.
            green_list = ["hasCalculatedProperty", "wasCalculatedBy", "hasValue"]
            ssplit = string3.split(".")
            if len(ssplit) == 2:
                if (ssplit[0] == "cmso") and (ssplit[1] not in green_list):
                    plot = False
            if string3 == "subClassOf":
                plot = False
            ssplit = string2.split(".")
            if string3 == "type":
                if (ssplit[0] == "cmso") and (ssplit[1] not in ["CalculatedProperty"]):
                    plot = False
                if (ssplit[0] == "cmso") and (ssplit[1] in ["AtomicScaleSample"]):
                    dot.node(
                        string1,
                        label=string1,
                        shape=styledict[istype1]["shape"],
                        style=styledict[istype1]["style"],
                        color=styledict[istype1]["color"],
                        fontsize=styledict[istype1]["fontsize"],
                        fontname=styledict[istype1]["fontname"],
                    )
                    plot = False
        
        elif sample_view:
            green_list = ['wasDerivedFrom', 'wasGeneratedBy']
            if string3 not in green_list:
                plot = False
            

        if hide_types and (string3 == "type"):
            plot = False

        if not plot:
            continue

        if istype1 == "Literal":
            id1 = str(uuid.uuid4())
        else:
            id1 = string1
        dot.node(
            id1,
            label=string1,
            shape=styledict[istype1]["shape"],
            style=styledict[istype1]["style"],
            color=styledict[istype1]["color"],
            fontsize=styledict[istype1]["fontsize"],
            fontname=styledict[istype1]["fontname"],
        )

        if istype2 == "Literal":
            id2 = str(uuid.uuid4())
        else:
            id2 = string2
        dot.node(
            id2,
            label=string2,
            shape=styledict[istype2]["shape"],
            style=styledict[istype2]["style"],
            color=styledict[istype2]["color"],
            fontsize=styledict[istype2]["fontsize"],
            fontname=styledict[istype2]["fontname"],
        )

        dot.edge(
            id1,
            id2,
            color=styledict["edgecolor"],
            label=string3,
            fontsize=styledict[istype2]["fontsize"],
            fontname=styledict[istype2]["fontname"],
        )

    return dot

def _id(item):
    return str(item).replace(':', '_')

def visualize_provenance(
    prov,
    rankdir="TB",
    size=None,
    layout="dot",
):
    dot = graphviz.Digraph()
    dot.attr(
        rankdir=rankdir,
        style="filled",
        size=size,
        layout=layout,
        overlap="false",
    )
    #add all nodes
    for key in prov.keys():
        nid = _id(key)
        #if "activity" in key:
        dot.node(nid, label=prov[key]['label'], 
                shape='box', 
                color="#C9DAF8", 
                style="filled",
                fontname='Helvetica',
                fontsize='8')
        #else:
        #    dot.node(nid, label=prov[key]['label'], 
        #            shape='parallelogram', 
        #            color="#C9DAF8", 
        #            style="filled",
        #            fontname='Helvetica',
        #            fontsize='8')
    #add all edges
    for key, val in prov.items():
        if 'inputs' in val.keys():
            if val['operation'] == 'input_parameter':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='input_param', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            if val['operation'] == 'output_parameter':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='output_param', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            elif val['operation'] == 'sample_for_activity':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='input_sample', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            elif val['operation'] == 'sample_output':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='output_sample', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            else:
                operation_id = str(uuid.uuid4())
                operation = dot.node(operation_id, label=val['operation'], 
                                    color="#E6B8AF", 
                                    shape='box', 
                                    style='filled',
                                    fontname='Helvetica',
                                    fontsize='8')
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), operation_id, label='input', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
                dot.edge(operation_id, _id(key), label='output', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
    return dot