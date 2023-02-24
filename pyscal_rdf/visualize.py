import graphviz
import os

def get_title_from_BNode(x):
    return x.title()

def get_string_from_URI(x):
    raw = x.title().split("#")
    if len(raw)>1:
        return raw[-1]
    else:
        return ".".join(x.title().split("/")[-2:])

def parse_object(x):
    if isinstance(x, BNode):
        return get_title_from_BNode(x)
    elif isinstance(x, URIRef):
        return get_string_from_URI(x)
    
def visualize_graph(g, shape="box",
            style="filled", 
            node_color="#AB47BC",
            edge_color="#37474F"):
    
    dot = graphviz.Digraph()
    for k in g:
        string1 = parse_object(k[0])
        dot.node(string1, shape=shape, style=style, color=node_color)
        string2 = parse_object(k[2])
        dot.node(string2, shape=shape, style=style, color=node_color)
        string3 = parse_object(k[1])
        dot.edge(string2, string1, color=edge_color, label=string3)
    
    return dot