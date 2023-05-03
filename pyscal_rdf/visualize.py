import graphviz
import os
from rdflib import BNode, URIRef, Literal
import uuid
import json
import ipycytoscape

def get_title_from_BNode(x):
    return x.toPython()

def get_string_from_URI(x, ):
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
            backend="ipycytoscape",
            edge_color="#37474F",
            styledict=styledict,
            graph_attr ={'rankdir': 'LR'},
            layoutname='cola'):
    if backend=='ipycytoscape':
        return _visualize_graph_ipycytoscape_backend(g,
                                                edge_color=edge_color,
                                                styledict=styledict,
                                                layoutname=layoutname)
    else:
        return _visualize_graph_graphviz_backend(g,
                                                edge_color=edge_color,
                                                styledict=styledict,
                                                graph_attr=graph_attr)
    
def _visualize_graph_ipycytoscape_backend(g,
                                          edge_color="#37474F",
                                          styledict=styledict,
                                          layoutname='cola'):
    #first step is to create the json file
    # we can start with a dict
    gdict = {}
    gdict["nodes"] = []
    gdict["edges"] = []
    for k in g:
        string1, istype1 = parse_object(k[0])
        id1 = _fix_id(string1, istype1)
        d = {"data": {"id": str(id1), 
                      "label": str(string1), 
                      "classes": str(istype1),
                      "shape": _switch_box(styledict[istype1]['shape']),
                      "width": len(string1)*7,
                      "fontsize": styledict[istype1]['fontsize']}}
        gdict["nodes"].append(d)

        string2, istype2 = parse_object(k[2])
        id2 = _fix_id(string2, istype2)
        d = {"data": {"id": str(id2), 
                      "label": str(string2), 
                      "classes": str(istype2),
                      "shape": _switch_box(styledict[istype2]['shape']),
                      "width": len(string2)*7,
                      "fontsize": styledict[istype1]['fontsize']}}
        gdict["nodes"].append(d)
        
        string3, istype3 = parse_object(k[1])
        id3 = str(uuid.uuid4())
        d = {"data": {"id": str(id3), 
                      "label": str(string3), 
                      "source": id1, 
                      "target": id2,
                      "fontsize": styledict[istype1]['fontsize']}}
        gdict["edges"].append(d)
    
    #graph is complete
    #lets try without any customisation
    my_style = [
        {'selector': 'node',
         'style': {
            'font-family': 'arial',
            'font-size': '10px',
            'label': 'data(label)',
            'shape': 'data(shape)',
            "text-valign": "center",
            "text-halign": "center",
            "width": 'data(width)',
            "font-size": 'data(fontsize)',
         }}, 
        {'selector': 
            'edge',
            'style': {
                'font-family': 'arial',
                'label': 'data(label)',
                "font-size": 'data(fontsize)',
                "target-arrow-shape": "triangle",
                "target-arrow-color": edge_color,
                "curve-style": "bezier",
            }},
         {'selector': 
            'node[classes="BNode"]',
            'style': {
                'background-color': styledict['BNode']['color'],
            }}, 
         {'selector': 
            'node[classes="URIRef"]',
            'style': {
                'background-color': styledict['URIRef']['color'],
            }}, 
         {'selector': 
            'node[classes="Literal"]',
            'style': {
                'background-color': styledict['Literal']['color']
            }}, 
    ]
      
    ipycyobj = ipycytoscape.CytoscapeWidget()
    ipycyobj.set_layout(
        name=layoutname,
        avoidOverlap=True,
        animate=True,)
    ipycyobj.graph.add_graph_from_json(gdict, directed=True)
    ipycyobj.set_style(my_style)
    return ipycyobj   
    

def _visualize_graph_graphviz_backend(g,  
            edge_color="#37474F",
            styledict=styledict,
            graph_attr ={'rankdir': 'LR'},
            rankdir='LR'):
    
    dot = graphviz.Digraph()
    for key, val in graph_attr.items():
        dot.graph_attr[key] = val
    for k in g:
        string1, istype1 = parse_object(k[0])
        if istype1 == 'Literal':
            id1 = str(uuid.uuid4())
        else:
            id1 = string1 
        dot.node(id1, label=string1, shape=styledict[istype1]["shape"], 
                 style=styledict[istype1]["style"], 
                 color=styledict[istype1]["color"],
                 fontsize=styledict[istype1]["fontsize"])
        
        string2, istype2 = parse_object(k[2])
        if istype2 == 'Literal':
            id2 = str(uuid.uuid4())
        else:
            id2 = string2 
        dot.node(id2, label=string2, shape=styledict[istype2]["shape"], 
                 style=styledict[istype2]["style"], 
                 color=styledict[istype2]["color"],
                 fontsize=styledict[istype2]["fontsize"])
        
        string3, istype = parse_object(k[1])
        dot.edge(id1, id2, color=edge_color, label=string3, fontsize=styledict[istype2]["fontsize"])
    
    return dot