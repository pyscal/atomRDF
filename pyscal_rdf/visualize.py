import graphviz
import os
from rdflib import BNode, URIRef, Literal, Namespace
import uuid
import json
import ipycytoscape


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
            backend="ipycytoscape",
            edge_color="#37474F",
            styledict=styledict,
            graph_attr ={'rankdir': 'LR'},
            layoutname='cola',
            hide_types=False,
            workflow_view=False):
    if backend=='ipycytoscape':
        return _visualize_graph_ipycytoscape_backend(g,
                                                edge_color=edge_color,
                                                styledict=styledict,
                                                layoutname=layoutname,
                                                hide_types=hide_types,
                                                workflow_view=workflow_view)
    else:
        return _visualize_graph_graphviz_backend(g,
                                                edge_color=edge_color,
                                                styledict=styledict,
                                                graph_attr=graph_attr,
                                                hide_types=hide_types,
                                                workflow_view=workflow_view)
    
def _visualize_graph_ipycytoscape_backend(g,
                                          edge_color="#37474F",
                                          styledict=styledict,
                                          layoutname='cola',
                                          hide_types=False,
                                          workflow_view=False):
    #first step is to create the json file
    # we can start with a dict
    gdict = {}
    gdict["nodes"] = []
    gdict["edges"] = []
    for k in g:
        string1, istype1 = parse_object(k[0])
        string2, istype2 = parse_object(k[2])
        string3, istype3 = parse_object(k[1])

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
                             fontsize=styledict[istype1]["fontsize"])
                    plot=False

        if hide_types and (string3 == 'type'):
            plot = False

        if not plot:
            continue

        id1 = _fix_id(string1, istype1)
        d = {"data": {"id": str(id1), 
                      "label": str(string1), 
                      "classes": str(istype1),
                      "shape": _switch_box(styledict[istype1]['shape']),
                      "width": len(string1)*7,
                      "fontsize": styledict[istype1]['fontsize']}}
        gdict["nodes"].append(d)

        id2 = _fix_id(string2, istype2)
        d = {"data": {"id": str(id2), 
                      "label": str(string2), 
                      "classes": str(istype2),
                      "shape": _switch_box(styledict[istype2]['shape']),
                      "width": len(string2)*7,
                      "fontsize": styledict[istype1]['fontsize']}}
        gdict["nodes"].append(d)
                
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
            rankdir='LR',
            hide_types=False,
            workflow_view=False):
    
    dot = graphviz.Digraph()
    for key, val in graph_attr.items():
        dot.graph_attr[key] = val
    
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
                             fontsize=styledict[istype1]["fontsize"])
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
                 fontsize=styledict[istype1]["fontsize"])
        
        if istype2 == 'Literal':
            id2 = str(uuid.uuid4())
        else:
            id2 = string2 
        dot.node(id2, label=string2, shape=styledict[istype2]["shape"], 
                 style=styledict[istype2]["style"], 
                 color=styledict[istype2]["color"],
                 fontsize=styledict[istype2]["fontsize"])
        
        dot.edge(id1, id2, color=edge_color, label=string3, fontsize=styledict[istype2]["fontsize"])
    
    return dot