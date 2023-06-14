from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header
from pyscal.crystal_structures import elements, structures
st_list = list(structures.keys())
el_list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Ar", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Ni", "Co", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "I", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

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
        self.element_box = dropdown('Element', 
                                 el_list, 
                                 value="Fe", 
                                 theme=theme)
        self.element_box_2 = dropdown('Element', 
                                 el_list, 
                                 value="Al", 
                                 theme=theme,
                                 disabled=True)
        
        self.text = "Create a bulk structure"
        self.header = header(self.text, theme=theme)
        self.panel = VBox(children=[self.header,
                                   self.dropdown,
                                   self.checkbox,
                                   self.repetition_box,
                                   self.lattice_parameter_box,
                                   self.element_box,
                                   self.element_box_2, 
                                   self.run_button])
    