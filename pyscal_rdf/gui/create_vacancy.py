from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header
from pyscal.crystal_structures import elements, structures
st_list = list(structures.keys())

class CreateVacancy:
    def __init__(self, theme="teal"):
        #create widgets
        self.vacancy_dropdown = dropdown('Add vacancy', 
                                ["concentration", "number"], 
                                theme=theme,
                                )
        self.vacancy_box = textbox('', 
                                value="0", 
                                theme=theme,
                                dtype="text")
        self.run_button = button('Generate', 
                                 theme=theme,
                                 disabled=True,
                                 tooltip="Add either vacancy concentration or number")
        self.panel = VBox(children=[HBox(children=[self.vacancy_dropdown, self.vacancy_box]), self.run_button])
