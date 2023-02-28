from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

from pyscal_rdf.visualize import visualize_graph

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")

styledict = {
    "BNode": {"color": "#ffe6ff", "shape": "box", "style": "filled"},
    "URIRef": {"color": "#ffffcc", "shape": "box", "style": "filled"},
    "Literal": {"color": "#e6ffcc", "shape": "parallelogram", "style": "filled"},
}

class StructureGraph(Graph):
    def __init__(self, structure_dict):
        super().__init__()
        self.sysdict = structure_dict
        self.sample = None
        self.material = None
    
    def create_graph(self, names=False, name_index="01"):
        if names:
            name_list = [f'{name_index}_Sample', f'{name_index}_Material',
                        f'{name_index}_ChemicalComposition', f'{name_index}_SimulationCell',
                        f'{name_index}_SimulationCell', f'{name_index}_CrystalStructure',
                        f'{name_index}_SpaceGroup', f'{name_index}_UnitCell',
                        f'{name_index}_UnitCell',]
        else:
            name_list = [None, None,
                        None, None,
                        None, None,
                        None, None,
                        None]
        self.add_sample(name=name_list[0])
        self.add_material(name=name_list[1])
        self.add_chemical_composition(name=name_list[2])
        self.add_simulation_cell(name=name_list[3])
        self.add_simulation_cell_properties(name=name_list[4])
        self.add_crystal_structure(name=name_list[5])
        self.add_space_group(name=name_list[6])
        self.add_unit_cell(name=name_list[7])
        self.add_lattice_properties(name=name_list[8])
        
    def add_sample(self, name=None):
        sample_01 = BNode(name)
        self.add((sample_01, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample_01
    
    def add_material(self, name=None):
        material_01 = BNode(name)
        self.add((self.sample, CMSO.hasMaterial, material_01))
        self.add((material_01, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material_01
    
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
        
    
    def add_simulation_cell_properties(self, name=None):
        if name is not None:
            name = f'{name}Length'
        simulation_cell_length_01 = BNode(name)
        self.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length_01))
        self.add((simulation_cell_length_01, RDF.type, CMSO.SimulationCellLength))
        self.add((simulation_cell_length_01, CMSO.hasLength_x, Literal(self.sysdict["SimulationCellLengthX"], datatype=XSD.float)))
        self.add((simulation_cell_length_01, CMSO.hasLength_y, Literal(self.sysdict["SimulationCellLengthY"], datatype=XSD.float)))
        self.add((simulation_cell_length_01, CMSO.hasLength_z, Literal(self.sysdict["SimulationCellLengthZ"], datatype=XSD.float)))
        
        if name is not None:
            name = f'{name}Vector01'
        simulation_cell_vector_01 = BNode(name)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        self.add((simulation_cell_vector_01, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorA"][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorA"][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorA"][2], datatype=XSD.float)))

        if name is not None:
            name = f'{name}Vector02'
        simulation_cell_vector_02 = BNode(name)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        self.add((simulation_cell_vector_02, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorB"][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorB"][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorB"][2], datatype=XSD.float)))

        if name is not None:
            name = f'{name}Vector03'
        simulation_cell_vector_03 = BNode(name)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        self.add((simulation_cell_vector_03, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorC"][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorC"][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorC"][2], datatype=XSD.float)))

        if name is not None:
            name = f'{name}Angle'
        simulation_cell_angle_01 = BNode(name)
        self.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle_01))
        self.add((simulation_cell_angle_01, RDF.type, CMSO.SimulationCellAngle))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_alpha, Literal(self.sysdict["SimulationCellAngleAlpha"], datatype=XSD.float)))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_beta, Literal(self.sysdict["SimulationCellAngleBeta"], datatype=XSD.float)))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_gamma, Literal(self.sysdict["SimulationCellAngleGamma"], datatype=XSD.float)))
        
    
    def add_crystal_structure(self, name=None):
        crystal_structure_01 = BNode(name)
        self.add((self.material, CMSO.hasStructure, crystal_structure_01))
        self.add((crystal_structure_01, RDF.type, CMSO.CrystalStructure))
        self.add((crystal_structure_01, CMSO.hasAltName, Literal(self.sysdict["CrystalStructureName"], datatype=XSD.string)))
        self.crystal_structure = crystal_structure_01
        
    def add_space_group(self, name=None):
        space_group_01 = BNode(name)
        self.add((self.material, CMSO.hasSpaceGroup, space_group_01))
        self.add((space_group_01, RDF.type, CMSO.SpaceGroup))
        self.add((space_group_01, CMSO.hasSpaceGroupSymbol, Literal(self.sysdict["SpaceGroupSymbol"], datatype=XSD.string)))
        self.add((space_group_01, CMSO.hasSpaceGroupNumber, Literal(self.sysdict["SpaceGroupNumber"], datatype=XSD.integer)))
    
    def add_unit_cell(self, name=None):
        unit_cell_01 = BNode(name)
        self.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell_01))
        self.add((unit_cell_01, RDF.type, CMSO.UnitCell))
        self.unit_cell = unit_cell_01
        
    def add_lattice_properties(self, name=None):
        if name is not None:
            name = f'{name}LatticeParameter'
        lattice_parameter_01 = BNode(name)
        self.add((self.unit_cell, CMSO.hasLatticeParamter, lattice_parameter_01))
        self.add((lattice_parameter_01, RDF.type, CMSO.LatticeParameter))
        self.add((lattice_parameter_01, CMSO.hasLength_x, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))
        self.add((lattice_parameter_01, CMSO.hasLength_y, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))
        self.add((lattice_parameter_01, CMSO.hasLength_z, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))

        if name is not None:
            name = f'{name}LatticeAngle'
        lattice_angle_01 = BNode(name)
        self.add((self.unit_cell, CMSO.hasAngle, lattice_angle_01))
        self.add((lattice_angle_01, RDF.type, CMSO.LatticeAngle))
        self.add((lattice_angle_01, CMSO.hasAngle_alpha, Literal(90, datatype=XSD.float)))
        self.add((lattice_angle_01, CMSO.hasAngle_beta, Literal(90, datatype=XSD.float)))
        self.add((lattice_angle_01, CMSO.hasAngle_gamma, Literal(90, datatype=XSD.float)))        
        
        
    def visualise(self, edge_color="#37474F",
            styledict=styledict):
        return visualize_graph(self, edge_color=edge_color,
            styledict=styledict)
    
    def write(filename):
        with open(filename, "w") as fout:
            fout.write(self.serialize(format="json-ld"))        