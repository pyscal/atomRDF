from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header
from pyscal.crystal_structures import elements, structures
st_list = list(structures.keys())

class CreateStructure:
    def __init__(self, theme="teal"):
        #create widgets
        self.dropdown = dropdown('Structure', 
                                 st_list, 
                                 value="bcc", 
                                 theme=theme)
        self.run_button = button('Generate', 
                                 theme=theme)
        self.checkbox = checkbox('Add to graph', 
                                 theme=theme)
        self.repetition_box = textbox('Repetitions', 
                                      value=1, 
                                      theme=theme,
                                     dtype="int")
        self.lattice_parameter_box = textbox('LatticeParameter', 
                                      value=1, 
                                      theme=theme,
                                      dtype="float")
        self.element_box = textbox('Element', 
                                      value="Fe", 
                                      theme=theme,
                                      dtype="text")
        self.element_box_2 = textbox('Element', 
                                      value="Al", 
                                      theme=theme,
                                      dtype="text",
                                      disabled=True)
        self.text = "Create a lattice"
        self.header = header(self.text, theme=theme)
        self.panel = VBox(children=[self.header,
                                   self.dropdown,
                                   self.checkbox,
                                   self.repetition_box,
                                   self.lattice_parameter_box,
                                   self.element_box,
                                   self.element_box_2, 
                                   self.run_button])
    