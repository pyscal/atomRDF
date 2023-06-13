from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header
from pyscal.crystal_structures import elements, structures
st_list = list(structures.keys())

class CreateGB:
    def __init__(self, theme="teal"):
        #create widgets
        self.axis_box = textbox('Axis', 
                                value="0 0 1", 
                                theme=theme,
                                dtype="text",)
        self.gb_plane_box = textbox('GB Plane', 
                                value="3 -1 0", 
                                theme=theme,
                                dtype="text")
        self.rep_box = textbox('Repetitions', 
                                value="1 1 1", 
                                theme=theme,
                                dtype="text")
        #self.sigma_dropdown = dropdown('Sigma', 
        #                              value=5, 
        #                              theme=theme,
        #                             dtype="int")
        self.sigma_button = button('Find sigma', 
                                 theme=theme,
                                 disabled=False,
                                 tooltip="Calculate possible sigma values, only 10 results are shown")
        self.gb_plane_button = button('Find possible GB planes', 
                                 theme=theme,
                                 disabled=False,
                                 tooltip="Calculate possible gb plane, only 10 results are shown")
        self.run_button = button('Generate', 
                                 theme=theme,
                                 disabled=True,
                                 tooltip="Creating a structure using the above methods is necessary. Reading in a file does not work")
        self.widget_list = [self.axis_box, self.sigma_button]
        self.panel = HBox(children=self.widget_list)
        self.sigma_dropdown = None
        self.sigma_hb = None
        self.gb_plane_dropdown = None
        self.gb_plane_hb = None
