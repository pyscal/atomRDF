import graphviz
import os
from rdflib import BNode, URIRef, Literal, Namespace
import uuid
import json


def get_title_from_BNode(x):
    return x.toPython()

def get_string_from_URI(x):
    """
    Extract a presentable string from URI

    Also differentiate between fixed notes and URIs, and assign color
    """
    raw = x.toPython()
    #first try splitting by #
    rawsplit = raw.split("#")
    if len(rawsplit) > 1:
        return rawsplit[-1], "URIRef"
    
    #try splitting by = for chebi values
    if 'CHEBI' in raw:
        rawsplit = raw.split("=")
        rawsplit = rawsplit[-1].split(":")
        if len(rawsplit) > 1:
            return ".".join(rawsplit[-2:]), "URIRef"
    
    if 'sample:' in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    #just a normal url split now
    rawsplit = raw.split("/")
    if len(rawsplit) > 1:
        return ".".join(rawsplit[-2:]), "URIRef"

    rawsplit = raw.split(':')
    if len(rawsplit) == 2:
        return "_".join(rawsplit), "BNode"        

    #none of the conditions, worked, which means its a hex string
    return raw, "BNode"

def parse_object(x):
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
        return 'rectangle'
    #remember that only boxes will be used, circles no!
    
def _fix_id(string1, istype1):
    if istype1 == 'Literal':
        id1 = str(uuid.uuid4())
    else:
        id1 = string1
    return id1

def visualize_graph(g,
            styledict=styledict,
            rankdir='TB',
            hide_types=False,
            workflow_view=False,
            size=None,
            layout='dot'):
    
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
            #we collapse sample information
            #if cmso.connector is found, only use it is it is cmso.hasCalculated
            #all sub sample props, indicated by sample_x_jsjsj will be ignored.
            ssplit = string3.split('.')
            if (len(ssplit) == 2):
                if (ssplit[0] == 'cmso') and (ssplit[1] != "hasCalculatedProperty"):
                    plot = False
            if string3 == 'subClassOf':
                plot = False
            ssplit = string2.split('.')
            if string3 == 'type':
                if (ssplit[0] == 'cmso') and (ssplit[1] not in  ["CalculatedProperty"]):
                    plot = False
                if (ssplit[0] == 'cmso') and (ssplit[1] in  ["AtomicScaleSample"]):
                    dot.node(string1, label=string1, shape=styledict[istype1]["shape"], 
                             style=styledict[istype1]["style"], 
                             color=styledict[istype1]["color"],
                             fontsize=styledict[istype1]["fontsize"],
                             fontname=styledict[istype1]["fontname"])
                    plot=False

        if hide_types and (string3 == 'type'):
            plot = False

        if not plot:
            continue

        if istype1 == 'Literal':
            id1 = str(uuid.uuid4())
        else:
            id1 = string1 
        dot.node(id1, label=string1, shape=styledict[istype1]["shape"], 
                 style=styledict[istype1]["style"], 
                 color=styledict[istype1]["color"],
                 fontsize=styledict[istype1]["fontsize"],
                 fontname=styledict[istype1]["fontname"])
        
        if istype2 == 'Literal':
            id2 = str(uuid.uuid4())
        else:
            id2 = string2 
        dot.node(id2, label=string2, shape=styledict[istype2]["shape"], 
                 style=styledict[istype2]["style"], 
                 color=styledict[istype2]["color"],
                 fontsize=styledict[istype2]["fontsize"],
                 fontname=styledict[istype2]["fontname"])
        
        dot.edge(id1, id2, color=styledict["edgecolor"], 
            label=string3, 
            fontsize=styledict[istype2]["fontsize"],
            fontname=styledict[istype2]["fontname"])
    
    return dot