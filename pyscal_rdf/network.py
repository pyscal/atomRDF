import networkx as nx
import matplotlib.pyplot as plt


class Network:
    """
    Network representation of CMSO
    """
    def __init__(self):
        self.g = nx.DiGraph()
    
    def add(self, sub, pred, obj, dtype=None):
        self.g.add_node(sub)
        self.g.add_node(pred)
        self.g.add_node(obj, dtype=dtype)            
        self.g.add_edge(sub, pred)
        self.g.add_edge(pred, obj)
    
    def draw(self):
        nx.draw(self.g, with_labels=True, font_weight='bold')
        
    def get_shortest_path(self, source, destination):
        path = nx.shortest_path(self.g, source=source, target=destination)
        return path
        