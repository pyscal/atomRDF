from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header, upload


class WriteStructure:
    def __init__(self, theme="teal"):
        #create widgets
        self.write_formats = ['lammps-data', 
                             'lammps-dump',
                              'poscar',
                             'turtle',
                             'json-ld',
                             'xml',
                             'n3']
        self.dropdown = dropdown('File format',
                                self.write_formats,
                                value="lammps-data")
        self.run_button = button('Generate', 
                                 theme=theme)
        self.text = "Download file"
        self.header = header(self.text, theme=theme)
        self.output = output()
        self.panel = VBox(children=[self.header,
                                   HBox(children=[self.dropdown,
                                                 self.run_button]),
                                   self.output])

