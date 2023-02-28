import graphviz
import os
from rdflib import BNode, URIRef, Literal

def get_title_from_BNode(x):
    return x.toPython()

def get_string_from_URI(x):
    raw = x.toPython().split("#")
    if len(raw)>1:
        return raw[-1]
    else:
        return ".".join(x.toPython().split("/")[-2:])

def parse_object(x):
    if isinstance(x, BNode):
        return get_title_from_BNode(x), "BNode"
    elif isinstance(x, URIRef):
        return get_string_from_URI(x), "URIRef"
    elif isinstance(x, Literal):
        return str(x.title()), "Literal"

styledict = {
    "BNode": {"color": "#ffe6ff", "shape": "box", "style": "filled"},
    "URIRef": {"color": "#ffffcc", "shape": "box", "style": "filled"},
    "Literal": {"color": "#e6ffcc", "shape": "parallelogram", "style": "filled"},
}

def visualize_graph(g,  
            edge_color="#37474F",
            styledict=styledict):
    
    dot = graphviz.Digraph()
    for k in g:
        string1, istype1 = parse_object(k[0])
        dot.node(string1, shape=styledict[istype1]["shape"], 
                 style=styledict[istype1]["style"], 
                 color=styledict[istype1]["color"])
        
        string2, istype2 = parse_object(k[2])
        dot.node(string2, shape=styledict[istype2]["shape"], 
                 style=styledict[istype2]["style"], 
                 color=styledict[istype2]["color"])
        
        string3, istype = parse_object(k[1])
        dot.edge(string2, string1, color=edge_color, label=string3)
    
    return dot