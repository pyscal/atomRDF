from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

from pyscal_rdf.visualize import visualize_graph

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")

styledict = {
    "BNode": {"color": "#ffe6ff", "shape": "box", "style": "filled"},
    "URIRef": {"color": "#ffffcc", "shape": "box", "style": "filled"},
    "Literal": {"color": "#e6ffcc", "shape": "parallelogram", "style": "filled"},
}

class StructureGraph(Graph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sysdict = kwargs['structure_dict']
        self.sample = None
        self.material = None
    
    def add_sample(self, name=None):
        sample_01 = BNode(name)
        self.add((sample_01, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample_01
    
    def add_material(self, name=None):
        material_01 = BNode(name)
        self.add((sample_01, CMSO.hasMaterial, material_01))
        self.add((material_01, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material
    
    def add_chemical_composition(self, name=None):
        chem_comp = ["=".join([x, str(y)]) for x,y in zip(self.sysdict["ChemicalCompositionElement"], self.sysdict["ChemicalCompositionRatio"])]
        chemical_composition_01 = BNode(name)
        self.add((self.material, CMSO.hasComposition, chemical_composition_01))
        self.add((chemical_composition_01, RDF.type, CMSO.ChemicalComposition))
        self.add((chemical_composition_01, CMSO.hasElementRatio, Literal(chem_comp[0], datatype=XSD.string)))
        self.add((chemical_composition_01, CMSO.hasElementRatio, Literal(chem_comp[1], datatype=XSD.string)))
    
    def add_simulation_cell(self, name=None):
        simulation_cell_01 = BNode(name)
        self.add((self.sample, CMSO.hasSimulationCell, simulation_cell_01))
        self.add((simulation_cell_01, RDF.type, CMSO.SimulationCell))
        self.add((simulation_cell_01, CMSO.hasVolume, Literal(self.sysdict["CellVolume"], datatype=XSD.float)))
        self.add((self.sample, CMSO.hasNumberOfAtoms, Literal(self.sysdict["NumberOfAtoms"], datatype=XSD.integer)))
        self.simulation_cell = simulation_cell_01
    
    def add_simulation_cell_lengths(self, name=None)
    
    

        
        
    def visualise(self, edge_color="#37474F",
            styledict=styledict):
        return visualize_graph(self, edge_color=edge_color,
            styledict=styledict)