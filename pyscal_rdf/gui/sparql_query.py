from ipywidgets import HBox, VBox, Layout
from pyscal_rdf.gui.create import output, dropdown, button, checkbox, textbox, header, upload

class SparqlQuery:
    def __init__(self, theme="teal"):
        #create widgets
        self.query_box = textbox("", 
                                value="""
            PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
            SELECT DISTINCT ?symbol
            WHERE {
                ?sample cmso:hasNumberOfAtoms ?number .
                ?sample cmso:hasMaterial ?material .
                ?material cmso:hasStructure ?structure .
                ?structure cmso:hasSpaceGroup ?spacegroup .
                ?spacegroup cmso:hasSpaceGroupSymbol ?symbol .
            FILTER (?number="2"^^xsd:integer)
            }
            """,
            dtype="textarea")
        self.output = output(theme=theme)
        self.run_button = button('Query', 
                                 theme=theme)
        self.panel = VBox(children=[self.query_box,
                                   self.run_button,
                                   self.output])