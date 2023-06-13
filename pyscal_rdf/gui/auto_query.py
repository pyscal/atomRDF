from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header, upload, text
from pyscal_rdf.gui.properties import dataprops, options

class AutoQuery:
    def __init__(self, theme="teal"):
        self.first_line = text("Find all samples with..")
        self.property_dropdown = dropdown("", options, 
                                          value="Bravais lattice ")
        self.compare_options = dropdown("", [" equal to ", 
                                                   " in between ",
                                                   " equal to either "],
                                       value=" equal to ")
        self.compare_output = output(theme=theme)
        self.input_output = output(theme=theme)
        self.input_field_1 = textbox("", value="bcc", dtype="text")
        self.input_field_2 = textbox("", value="", dtype="text")
        self.run_button = button("Run query", theme=theme)
        self.show_button = button("Show query", theme=theme)
        self.result_output = output(theme=theme)
        self.plot_structure_button = button("Plot structure", theme=theme)
        self.visualise_graph_button = button("Visualise graph", theme=theme)
        self.result_plotter = output(theme=theme)
        self.w1 = HBox(children=[ self.property_dropdown, 
                       self.compare_output,
                       self.input_output])
        self.panel = VBox(children=[self.first_line, 
               self.w1,
               HBox(children=[self.run_button, self.show_button]),
               self.result_output,
               self.result_plotter])