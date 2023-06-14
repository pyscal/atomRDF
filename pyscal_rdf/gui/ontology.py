from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header, upload

class Ontology:
    def __init__(self, theme="teal"):
        #create widgets
        self.text = "pyscal-rdf uses the Computational Material Sample Ontology"
        self.header = header(self.text, theme=theme)
        self.output = output()
        self.panel = VBox(children=[self.header, 
                                    self.output])
