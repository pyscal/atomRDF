from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import dropdown, button, checkbox, textbox, header
from pyscal.crystal_structures import elements, structures
el_list = list(elements.keys())

class CreateElement:
    def __init__(self, theme="teal"):
        #create widgets
        self.dropdown = dropdown('Element', 
                                 el_list, 
                                 value="Fe", 
                                 theme=theme)
        self.run_button = button('Generate', 
                                 theme=theme)
        self.checkbox = checkbox('Add to graph', 
                                 theme=theme)
        self.repetition_box = textbox('Repetitions', 
                                      value=1, 
                                      theme=theme,
                                      dtype="int")
        self.text = "Create using element symbol"
        self.header = header(self.text, 
                             theme=theme, tooltip="create an element")
        self.panel = VBox(children=[self.header,
                                   self.dropdown,
                                   self.checkbox,
                                   self.repetition_box,
                                   self.run_button])
    