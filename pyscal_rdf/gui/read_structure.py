from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import dropdown, button, checkbox, textbox, header, upload

class ReadStructure:
    def __init__(self, theme="teal"):
        #create widgets
        self.read_formats = ['lammps-data', 
                             'lammps-dump', 
                             'poscar', 
                             'cif']
        self.dropdown = dropdown('File format',
                                self.read_formats,
                                value="poscar")
        self.upload = upload()
        self.run_button = button('Generate', 
                                 theme=theme)
        self.text = "Read structure"
        self.header = header(self.text, theme=theme)
        self.panel = VBox(children=[self.header,
                                   HBox(children=[self.upload,
                                                 self.dropdown]),
                                   self.run_button])
        
        
