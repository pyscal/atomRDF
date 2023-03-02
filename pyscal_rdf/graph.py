from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS
from rdflib.store import NO_STORE, VALID_STORE

import os
import numpy as np

from pyscal_rdf.visualize import visualize_graph
from pyscal_rdf.rdfutils import convert_to_dict
from pyscal.core import System
from pyscal.atoms import Atoms
from pyscal.core import System

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")

styledict = {
    "BNode": {"color": "#ffe6ff", 
              "shape": "box", 
              "style": "filled",
              "fontsize": "8"},
    "URIRef": {"color": "#ffffcc", 
               "shape": "box", 
               "style": "filled",
               "fontsize": "8"},
    "Literal": {"color": "#e6ffcc", 
                "shape": "parallelogram", 
                "style": "filled",
                "fontsize": "8"},
}

class StructureGraph:
    def __init__(self, graph_file=None):
        self.graph = Graph()
        self.graph.bind("cmso", CMSO)
        if graph_file is not None:
            if os.path.exists(graph_file):
                self.graph.parse(input_info)
        self.sample = None
        self.material = None
    
    def process_structure(self, structure):
        if isinstance(structure, System):
            self.sysdict = convert_to_dict(structure)
        elif os.path.exists(structure):
            sys = System(structure, format=format)
            self.sysdict = convert_to_dict(sys)
        
    def add_structure_to_graph(self, structure, names=False, name_index="01", format=None):
        self.process_structure(structure)
        #now add to graph
        self.create_graph(names=names, name_index=name_index)
        
    
    def create_graph(self, names=False, name_index="01"):
        if names:
            name_list = [f'{name_index}_Sample', f'{name_index}_Material',
                        f'{name_index}_ChemicalComposition', f'{name_index}_SimulationCell',
                        f'{name_index}_SimulationCell', f'{name_index}_CrystalStructure',
                        f'{name_index}_SpaceGroup', f'{name_index}_UnitCell',
                        f'{name_index}_UnitCell']
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
        self.add_atoms()
        
    def add_sample(self, name=None):
        sample_01 = BNode(name)
        self.graph.add((sample_01, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample_01
    
    def add_material(self, name=None):
        material_01 = BNode(name)
        self.graph.add((self.sample, CMSO.hasMaterial, material_01))
        self.graph.add((material_01, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material_01
    
    def add_chemical_composition(self, name=None):
        chem_comp = ["=".join([x, str(y)]) for x,y in zip(self.sysdict["ChemicalCompositionElement"], self.sysdict["ChemicalCompositionRatio"])]
        chemical_composition_01 = BNode(name)
        self.graph.add((self.material, CMSO.hasComposition, chemical_composition_01))
        self.graph.add((chemical_composition_01, RDF.type, CMSO.ChemicalComposition))
        for x in range(len(chem_comp)):
            self.graph.add((chemical_composition_01, CMSO.hasElementRatio, Literal(chem_comp[x], datatype=XSD.string)))
    
    def add_simulation_cell(self, name=None):
        simulation_cell_01 = BNode(name)
        self.graph.add((self.sample, CMSO.hasSimulationCell, simulation_cell_01))
        self.graph.add((simulation_cell_01, RDF.type, CMSO.SimulationCell))
        self.graph.add((simulation_cell_01, CMSO.hasVolume, Literal(np.round(self.sysdict["CellVolume"], decimals=2), datatype=XSD.float)))
        self.graph.add((self.sample, CMSO.hasNumberOfAtoms, Literal(self.sysdict["NumberOfAtoms"], datatype=XSD.integer)))
        self.simulation_cell = simulation_cell_01
        
    
    def add_simulation_cell_properties(self, name=None):
        uname = None
        if name is not None:
            uname = f'{name}Length'
        simulation_cell_length_01 = BNode(uname)
        self.graph.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length_01))
        self.graph.add((simulation_cell_length_01, RDF.type, CMSO.SimulationCellLength))
        self.graph.add((simulation_cell_length_01, CMSO.hasLength_x, Literal(self.sysdict["SimulationCellLengthX"], datatype=XSD.float)))
        self.graph.add((simulation_cell_length_01, CMSO.hasLength_y, Literal(self.sysdict["SimulationCellLengthY"], datatype=XSD.float)))
        self.graph.add((simulation_cell_length_01, CMSO.hasLength_z, Literal(self.sysdict["SimulationCellLengthZ"], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector01'
        simulation_cell_vector_01 = BNode(uname)
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        self.graph.add((simulation_cell_vector_01, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorA"][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorA"][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorA"][2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector02'
        simulation_cell_vector_02 = BNode(uname)
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        self.graph.add((simulation_cell_vector_02, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorB"][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorB"][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorB"][2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector03'
        simulation_cell_vector_03 = BNode(uname)
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        self.graph.add((simulation_cell_vector_03, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_x, Literal(self.sysdict["SimulationCellVectorC"][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_y, Literal(self.sysdict["SimulationCellVectorC"][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_z, Literal(self.sysdict["SimulationCellVectorC"][2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Angle'
        simulation_cell_angle_01 = BNode(uname)
        self.graph.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle_01))
        self.graph.add((simulation_cell_angle_01, RDF.type, CMSO.SimulationCellAngle))
        self.graph.add((simulation_cell_angle_01, CMSO.hasAngle_alpha, Literal(self.sysdict["SimulationCellAngleAlpha"], datatype=XSD.float)))
        self.graph.add((simulation_cell_angle_01, CMSO.hasAngle_beta, Literal(self.sysdict["SimulationCellAngleBeta"], datatype=XSD.float)))
        self.graph.add((simulation_cell_angle_01, CMSO.hasAngle_gamma, Literal(self.sysdict["SimulationCellAngleGamma"], datatype=XSD.float)))
        
    
    def add_crystal_structure(self, name=None):
        crystal_structure_01 = BNode(name)
        self.graph.add((self.material, CMSO.hasStructure, crystal_structure_01))
        self.graph.add((crystal_structure_01, RDF.type, CMSO.CrystalStructure))
        self.graph.add((crystal_structure_01, CMSO.hasAltName, Literal(self.sysdict["CrystalStructureName"], datatype=XSD.string)))
        self.crystal_structure = crystal_structure_01
        
    def add_space_group(self, name=None):
        space_group_01 = BNode(name)
        self.graph.add((self.crystal_structure, CMSO.hasSpaceGroup, space_group_01))
        self.graph.add((space_group_01, RDF.type, CMSO.SpaceGroup))
        self.graph.add((space_group_01, CMSO.hasSpaceGroupSymbol, Literal(self.sysdict["SpaceGroupSymbol"], datatype=XSD.string)))
        self.graph.add((space_group_01, CMSO.hasSpaceGroupNumber, Literal(self.sysdict["SpaceGroupNumber"], datatype=XSD.integer)))
    
            
    def add_unit_cell(self, name=None):
        unit_cell_01 = BNode(name)
        self.graph.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell_01))
        self.graph.add((unit_cell_01, RDF.type, CMSO.UnitCell))
        self.unit_cell = unit_cell_01
        
        #add bravais lattice
        uname = None
        if name is not None:
            uname = f'Bravais{name}'
        bravaislattice = BNode(uname)
        self.graph.add((self.unit_cell, CMSO.hasLattice, bravaislattice))
        self.graph.add((bravaislattice, RDF.type, CMSO.BravaisLattice))
        self.graph.add((bravaislattice, CMSO.hasLatticeSystem, Literal(self.sysdict["BravaisLattice"], datatype=XSD.string)))
        
    def add_lattice_properties(self, name=None):
        uname = None
        if name is not None:
            uname = f'{name}LatticeParameter'
        lattice_parameter_01 = BNode(uname)
        self.graph.add((self.unit_cell, CMSO.hasLatticeParamter, lattice_parameter_01))
        self.graph.add((lattice_parameter_01, RDF.type, CMSO.LatticeParameter))
        self.graph.add((lattice_parameter_01, CMSO.hasLength_x, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))
        self.graph.add((lattice_parameter_01, CMSO.hasLength_y, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))
        self.graph.add((lattice_parameter_01, CMSO.hasLength_z, Literal(self.sysdict["LatticeParameter"], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}LatticeAngle'
        lattice_angle_01 = BNode(uname)
        self.graph.add((self.unit_cell, CMSO.hasAngle, lattice_angle_01))
        self.graph.add((lattice_angle_01, RDF.type, CMSO.LatticeAngle))
        self.graph.add((lattice_angle_01, CMSO.hasAngle_alpha, Literal(90, datatype=XSD.float)))
        self.graph.add((lattice_angle_01, CMSO.hasAngle_beta, Literal(90, datatype=XSD.float)))
        self.graph.add((lattice_angle_01, CMSO.hasAngle_gamma, Literal(90, datatype=XSD.float)))        
        
    def add_atoms(self, name=None):
        for x in range(len(self.sysdict["Positions"])):
            #create atom
            atom = BNode()
            self.graph.add((self.sample, CMSO.hasAtom, atom))
            self.graph.add((atom, RDF.type, CMSO.Atom))
            position = BNode()
            self.graph.add((atom, CMSO.hasPositionVector, position))
            self.graph.add((position, RDF.type, CMSO.PositionVector))
            self.graph.add((position, CMSO.hasComponent_x, Literal(self.sysdict["Positions"][x][0],
                                                                  datatype=XSD.float)))
            self.graph.add((position, CMSO.hasComponent_y, Literal(self.sysdict["Positions"][x][1],
                                                                  datatype=XSD.float)))
            self.graph.add((position, CMSO.hasComponent_z, Literal(self.sysdict["Positions"][x][2],
                                                                  datatype=XSD.float)))
            #now add coordination
            element = BNode()
            self.graph.add((atom, CMSO.hasElement, element))
            self.graph.add((element, RDF.type, CMSO.Element))
            self.graph.add((element, CMSO.hasSymbol, Literal(str(self.sysdict["Element"][x]),
                                                            datatype=XSD.string)))
            #finally occupancy
            self.graph.add((atom, CMSO.hasCoordinationNumber, Literal(self.sysdict["Coordination"][x],
                                                                     datatype=XSD.integer)))
        
    def visualise(self, edge_color="#37474F",
            styledict=styledict, rankdir='BT'):
        return visualize_graph(self.graph, edge_color=edge_color,
            styledict=styledict, rankdir=rankdir)
    
    def write(self, filename, format="json-ld"):
        with open(filename, "w") as fout:
            fout.write(self.graph.serialize(format=format))
            
    
    def get_system_from_sample(self, sample):
        simcell = self.graph.value(sample, CMSO.hasSimulationCell)
        cell_vectors = [[], [], []]

        for s in self.graph.triples((simcell, CMSO.hasVector, None)):
            cell_vectors[0].append(self.graph.value(s[2], CMSO.hasComponent_x).toPython())
            cell_vectors[1].append(self.graph.value(s[2], CMSO.hasComponent_y).toPython())
            cell_vectors[2].append(self.graph.value(s[2], CMSO.hasComponent_z).toPython())
        #cell_vectors

        positions = []
        species = []

        for atom in self.graph.triples((sample, CMSO.hasAtom, None)):
            vector = self.graph.value(atom[2], CMSO.hasPositionVector)
            pt = []
            pt.append(self.graph.value(vector, CMSO.hasComponent_x).toPython())
            pt.append(self.graph.value(vector, CMSO.hasComponent_y).toPython())
            pt.append(self.graph.value(vector, CMSO.hasComponent_z).toPython())
            element = self.graph.value(atom[2], CMSO.hasElement)
            species.append(self.graph.value(element, CMSO.hasSymbol).toPython())
            positions.append(pt)
        
        atoms = {"positions": positions, "species": species}
        at = Atoms()
        at.from_dict(atoms)
        sys = System()
        sys.box = cell_vectors
        sys.atoms = at
        
        return sys
    
    def to_file(self, sample, filename=None, format="lammps-dump"):
        if filename is None:
            filename = os.path.join(os.getcwd(), "out")
        sys = self.get_system_from_sample(sample)
        if format=="ase":
            return sys.to_ase()
        else:
            sys.to_file(filename, format=format)
            