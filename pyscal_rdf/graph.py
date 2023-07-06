"""
Graph module contains the basic RDFGraph object in pyscal_rdf. This object gets a structure
as an input and annotates it with the CMSO ontology (PLDO and PODO too as needed). The annotated
object is stored in triplets.
"""

from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS
from rdflib.store import NO_STORE, VALID_STORE

import os
import numpy as np
import inspect
from ase.io import write
import copy
import pandas as pd

from pyscal_rdf.visualize import visualize_graph
from pyscal_rdf.network import OntologyNetwork
from pyscal_rdf.rdfsystem import System
import pyscal_rdf.properties as prp
#from pyscal3.core import System
from pyscal3.atoms import Atoms

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("https://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")

defstyledict = {
    "BNode": {"color": "#ffe6ff", 
              "shape": "box", 
              "style": "filled",
              "fontsize": "8",
              "fontname": "Courier"},
    "URIRef": {"color": "#ffffcc", 
               "shape": "box", 
               "style": "filled",
               "fontsize": "8"},
    "Literal": {"color": "#e6ffcc", 
                "shape": "parallelogram", 
                "style": "filled",
                "fontsize": "8"},
}

def _replace_keys(refdict, indict):
    for key, val in indict.items():
        if key in refdict.keys():
            if isinstance(val, dict):
                _replace_keys(refdict[key], indict[key])
            else:
                refdict[key] = val
    return refdict

class RDFGraph:
    def __init__(self, graph_file=None, 
        store="Memory", 
        store_file=None,
        identifier="default_graph"):
        
        #owlfile = os.path.join(os.path.dirname(__file__), "data/cmso.owl")
        #self.graph.parse(owlfile, format='xml')
        if store == "Memory":
            self.graph = Graph(store="Memory", identifier=identifier)

        elif store=="SQLAlchemy":
            if store_file is None:
                raise ValueError("store file is needed if store is not memory")
            self.graph = Graph(store="SQLAlchemy", identifier=identifier)
            uri = Literal(f"sqlite:///{store_file}")
            self.graph.open(uri, create=True)

        elif inspect.isclass(type(store)):
            try:
                prpath = store.path
                dbfile = os.path.join(prpath, 'project.db')
                #now start sqlalchemy instance
                self.graph = Graph(store="SQLAlchemy", identifier=identifier)
                uri = Literal(f"sqlite:///{dbfile}")
                self.graph.open(uri, create=True)
            except:
                raise ValueError("store should be pyiron_project, SQLAlchemy, or Memory")
        
        else:
            raise ValueError("store should be pyiron_project, SQLAlchemy, or Memory")

        self.graph.bind("cmso", CMSO)
        self.graph.bind("pldo", PLDO)
        
        if graph_file is not None:
            if os.path.exists(graph_file):
                self.graph.parse(graph_file)
        
        self.sample = None
        self.material = None
        self.sysdict = None
        self.sgraph = None
        self._query_graph = OntologyNetwork()
        self._atom_ids = None
    
    def process_structure(self, structure, format=None):
        """
        Convert a given :py:class:`pyscal.core.System` to a data dictionary which can be used for annotation
        and storing the data in the RDF Graph.

        Parameters
        ----------
        structure: :py:class:`pyscal.core.System`
            input structure

        Returns
        -------
        None
        """
        if isinstance(structure, System):
            #self.sysdict = convert_to_dict(structure)
            self.sys = structure
        elif os.path.exists(structure):
            sys = System(structure, format=format)
            #self.sysdict = convert_to_dict(sys)
            self.sys = sys
    
    def data(self, key):
        #this method gets info directly from the dict
        if key=="ChemicalCompositionElement":
            return list(self.sys.composition.keys())
        elif key=="ChemicalCompositionRatio":
            return [val for key, val in self.sys.composition.items()]
        elif key=="CellVolume":
            return self.sys.volume
        elif key=="NumberOfAtoms":
            return self.sys.natoms
        elif key=="SimulationCellLengthX":
            return self.sys.box_dimensions[0]
        elif key=="SimulationCellLengthY":
            return self.sys.box_dimensions[1]
        elif key=="SimulationCellLengthZ":
            return self.sys.box_dimensions[2]

        elif key=="SimulationCellVectorA":
            return self.sys.box[0]
        elif key=="SimulationCellVectorB":
            return self.sys.box[1]
        elif key=="SimulationCellVectorC":
            return self.sys.box[2]
        
        elif key=="SimulationCellAngleAlpha":
            return prp.get_angle(self.sys.box[0], self.sys.box[1])
        elif key=="SimulationCellAngleBeta":
            return prp.get_angle(self.sys.box[1], self.sys.box[2])
        elif key=="SimulationCellAngleGamma": 
            return prp.get_angle(self.sys.box[2], self.sys.box[0])
        
        elif key=="Element":
            if self.sys.atoms.species[0] is not None:
                return self.sys.atoms.species
            else:
                return self.sys.atoms.types
        elif key=="Coordination":
            return prp.get_coordination(self.sys)
        elif key=="Positions":
            return self.sys.atoms.positions
        elif key=="LatticeParameter":
            return self.sys.atoms._lattice_constant      
        
        elif key=="SpaceGroupSymbol":
            if self.sys._structure_dict is not None:
                symbol, number = prp.get_space_group(self.sys)
                return symbol
        elif key=="SpaceGroupNumber":
            if self.sys._structure_dict is not None:
                symbol, number = prp.get_space_group(self.sys)
                return number
            else:
                return None

        elif key=="CrystalStructureName":
            if self.sys._structure_dict is not None:
                return self.sys.atoms._lattice 
            else:
                return None

        elif key=="BravaisLattice":
            if self.sys._structure_dict is not None:
                return prp.get_bravais_lattice(self.sys)
            else:
                return None

        elif key=="BasisPositions":
            if self.sys._structure_dict is not None:
                return self.sys._structure_dict['positions']
            else:
                return None

        elif key=="BasisOccupancy":
            if self.sys._structure_dict is not None:
                return prp.get_basis(self.sys)
            else:
                return None

        elif key=="LatticeVectors":
            if self.sys._structure_dict is not None:
                return prp.get_lattice_vector(self.sys)
            else:
                return None

        #if self.sysdict is not None:
        #    if key in self.sysdict:
        #        return self.sysdict[key]
        return None
    
    def add(self, triple):
        if str(triple[2].toPython()) != 'None':
            self.graph.add(triple)
        
    def add_structure_to_graph(self, structure, names=True, name_index=None, format=None):
        """
        Add a given :py:class:`pyscal.core.System` to the Graph object

        Parameters
        ----------
        structure: :py:class:`pyscal.core.System`
            input structure

        names: bool
            if True, alphanumeric names will be used instead of random BNodes

        Returns
        -------
        None
        """
        self.process_structure(structure, format=format)
        #now add to graph
        if name_index is None:
            name_index = self.n_samples + 1
        self.create_graph(names=names, name_index=name_index)
        structure.sample = self.sample
        structure._atom_ids = copy.copy(self._atom_ids)
        structure.graph = self
    
    def create_graph(self, names=False, name_index="1"):
        """
        Create the RDF Graph from the data stored

        Parameters
        ----------
        names: bool
            if True, alphanumeric names will be used instead of random BNodes

        name_index: string
            Prefix to be added to identifiers, default 01        
        
        Returns
        -------
        None
        """

        if names:
            name_list = [f'Sample_{name_index}', f'Material_{name_index}',
                        f'ChemicalComposition_{name_index}', f'SimulationCell_{name_index}',
                        f'SimulationCell_{name_index}', f'CrystalStructure_{name_index}',
                        f'SpaceGroup_{name_index}', f'UnitCell_{name_index}',
                        f'UnitCell_{name_index}', f'Atom_{name_index}']
        else:
            name_list = [None, None,
                        None, None,
                        None, None,
                        None, None,
                        None, None]
        self.add_sample(name=name_list[0])
        self.add_material(name=name_list[1])
        self.add_chemical_composition(name=name_list[2])
        self.add_simulation_cell(name=name_list[3])
        self.add_simulation_cell_properties(name=name_list[4])
        self.add_crystal_structure(name=name_list[5])
        self.add_space_group(name=name_list[6])
        self.add_unit_cell(name=name_list[7])
        self.add_lattice_properties(name=name_list[8])
        self.add_atoms(name=name_list[9])

        self.add((CMSO.SimulationCellLength, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.LatticeParameter, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.Length, CMSO.hasUnit, URIRef("https://qudt.org/2.1/vocab/unit#ANGSTROM")))
        
        self.add((CMSO.SimulationCellAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.LatticeAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.Angle, CMSO.hasUnit, URIRef("https://qudt.org/2.1/vocab/unit#DEG")))
        
        self.add((CMSO.LatticeVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.SimulationCellVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.PositionVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.Vector, CMSO.hasUnit, URIRef("https://qudt.org/2.1/vocab/unit#ANGSTROM")))
        
        
    def add_sample(self, name=None):
        """
        Add a CMSO Sample object

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        sample_01 = BNode(name)
        self.add((sample_01, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample_01
    
    def add_material(self, name=None):
        """
        Add a CMSO Material object

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        material_01 = BNode(name)
        self.add((self.sample, CMSO.hasMaterial, material_01))
        self.add((material_01, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material_01
    
    def add_chemical_composition(self, name=None):
        """
        Add chemical composition

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        chem_comp = ["=".join([str(x), str(y)]) for x,y in zip(self.data("ChemicalCompositionElement"), self.data("ChemicalCompositionRatio"))]
        for x in range(len(chem_comp)):
            self.add((self.material, CMSO.hasElementRatio, Literal(chem_comp[x], datatype=XSD.string)))
    
    def add_simulation_cell(self, name=None):
        """
        Add a CMSO SimulationCell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        simulation_cell_01 = BNode(name)
        self.add((self.sample, CMSO.hasSimulationCell, simulation_cell_01))
        self.add((simulation_cell_01, RDF.type, CMSO.SimulationCell))
        self.add((simulation_cell_01, CMSO.hasVolume, Literal(np.round(self.data("CellVolume"), decimals=2), datatype=XSD.float)))
        self.add((self.sample, CMSO.hasNumberOfAtoms, Literal(self.data("NumberOfAtoms"), datatype=XSD.integer)))
        self.simulation_cell = simulation_cell_01
        
    
    def add_simulation_cell_properties(self, name=None):
        """
        Add a CMSO SimulationCell properties such as SimulationCellLength,
        and Vectors.

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        uname = None
        if name is not None:
            uname = f'{name}Length'
        simulation_cell_length_01 = BNode(uname)
        self.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length_01))
        self.add((simulation_cell_length_01, RDF.type, CMSO.SimulationCellLength))
        self.add((simulation_cell_length_01, CMSO.hasLength_x, Literal(self.data("SimulationCellLengthX"), datatype=XSD.float)))
        self.add((simulation_cell_length_01, CMSO.hasLength_y, Literal(self.data("SimulationCellLengthY"), datatype=XSD.float)))
        self.add((simulation_cell_length_01, CMSO.hasLength_z, Literal(self.data("SimulationCellLengthZ"), datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector01'
        simulation_cell_vector_01 = BNode(uname)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        self.add((simulation_cell_vector_01, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_x, Literal(self.data("SimulationCellVectorA")[0], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_y, Literal(self.data("SimulationCellVectorA")[1], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_z, Literal(self.data("SimulationCellVectorA")[2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector02'
        simulation_cell_vector_02 = BNode(uname)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        self.add((simulation_cell_vector_02, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_x, Literal(self.data("SimulationCellVectorB")[0], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_y, Literal(self.data("SimulationCellVectorB")[1], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_z, Literal(self.data("SimulationCellVectorB")[2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Vector03'
        simulation_cell_vector_03 = BNode(uname)
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        self.add((simulation_cell_vector_03, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_x, Literal(self.data("SimulationCellVectorC")[0], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_y, Literal(self.data("SimulationCellVectorC")[1], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_z, Literal(self.data("SimulationCellVectorC")[2], datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}Angle'
        simulation_cell_angle_01 = BNode(uname)
        self.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle_01))
        self.add((simulation_cell_angle_01, RDF.type, CMSO.SimulationCellAngle))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_alpha, Literal(self.data("SimulationCellAngleAlpha"), datatype=XSD.float)))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_beta, Literal(self.data("SimulationCellAngleBeta"), datatype=XSD.float)))
        self.add((simulation_cell_angle_01, CMSO.hasAngle_gamma, Literal(self.data("SimulationCellAngleGamma"), datatype=XSD.float)))
        
    
    def add_crystal_structure(self, name=None):
        """
        Add a CMSO Crystal Structure

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        crystal_structure_01 = BNode(name)
        self.add((self.material, CMSO.hasStructure, crystal_structure_01))
        self.add((crystal_structure_01, RDF.type, CMSO.CrystalStructure))    
        self.add((crystal_structure_01, CMSO.hasAltName, Literal(self.data("CrystalStructureName"), datatype=XSD.string)))
        self.crystal_structure = crystal_structure_01
        
    def add_space_group(self, name=None):
        """
        Add a CMSO Space Group

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        self.add((self.crystal_structure, CMSO.hasSpaceGroupSymbol, Literal(self.data("SpaceGroupSymbol"), datatype=XSD.string)))
        self.add((self.crystal_structure, CMSO.hasSpaceGroupNumber, Literal(self.data("SpaceGroupNumber"), datatype=XSD.integer)))
    
            
    def add_unit_cell(self, name=None):
        """
        Add a CMSO Unit Cell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        unit_cell_01 = BNode(name)
        self.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell_01))
        self.add((unit_cell_01, RDF.type, CMSO.UnitCell))
        self.unit_cell = unit_cell_01
        
        #add bravais lattice
        bv = None
        if self.data("BravaisLattice") is not None:
            bv = URIRef(self.data("BravaisLattice"))
            self.add((self.unit_cell, CMSO.hasBravaisLattice, bv))
        
    def add_lattice_properties(self, name=None):
        """
        Add CMSO lattice properties such as Lattice Parameter,
        and its lengths and angles. 

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        uname = None
        if name is not None:
            uname = f'{name}LatticeParameter'
        lattice_parameter_01 = BNode(uname)
        self.add((self.unit_cell, CMSO.hasLatticeParamter, lattice_parameter_01))
        self.add((lattice_parameter_01, RDF.type, CMSO.LatticeParameter))
        self.add((lattice_parameter_01, CMSO.hasLength_x, Literal(self.data("LatticeParameter"), datatype=XSD.float)))
        self.add((lattice_parameter_01, CMSO.hasLength_y, Literal(self.data("LatticeParameter"), datatype=XSD.float)))
        self.add((lattice_parameter_01, CMSO.hasLength_z, Literal(self.data("LatticeParameter"), datatype=XSD.float)))
        
        uname = None
        if name is not None:
            uname = f'{name}LatticeAngle'
        lattice_angle_01 = BNode(uname)
        self.add((self.unit_cell, CMSO.hasAngle, lattice_angle_01))
        self.add((lattice_angle_01, RDF.type, CMSO.LatticeAngle))
        self.add((lattice_angle_01, CMSO.hasAngle_alpha, Literal(90, datatype=XSD.float)))
        self.add((lattice_angle_01, CMSO.hasAngle_beta, Literal(90, datatype=XSD.float)))
        self.add((lattice_angle_01, CMSO.hasAngle_gamma, Literal(90, datatype=XSD.float)))        
        
    def add_atoms(self, name=None):
        """
        Add Atoms including their species and positions

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        self._atom_ids = []

        positions = self.data("Positions")
        elements = self.data("Element")
        coordination = self.data("Coordination")
        for x in range(len(positions)):
            uname = None
            if name is not None:
                uname = f'{name}_{x}'            
            #create atom
            atom = BNode(uname)
            self._atom_ids.append(atom)
            self.add((self.sample, CMSO.hasAtom, atom))
            self.add((atom, RDF.type, CMSO.Atom))

            uname = None
            if name is not None:
                uname = f'{name}_{x}_Position'            
            position = BNode(uname)
            self.add((atom, CMSO.hasPositionVector, position))
            self.add((position, RDF.type, CMSO.PositionVector))
            self.add((position, CMSO.hasComponent_x, Literal(positions[x][0],
                                                                  datatype=XSD.float)))
            self.add((position, CMSO.hasComponent_y, Literal(positions[x][1],
                                                                  datatype=XSD.float)))
            self.add((position, CMSO.hasComponent_z, Literal(positions[x][2],
                                                                  datatype=XSD.float)))
            #now add coordination
            uname = None
            if name is not None:
                uname = f'{name}_{x}_Element'            
            element = BNode(uname)
            self.add((atom, CMSO.hasElement, element))
            self.add((element, RDF.type, CMSO.Element))
            self.add((element, CMSO.hasSymbol, Literal(str(elements[x]),
                                                            datatype=XSD.string)))
            #finally occupancy
            self.add((atom, CMSO.hasCoordinationNumber, Literal(coordination[x],
                                                                     datatype=XSD.integer)))

    
    
    def add_gb(self, gb_dict, name=None):
        """
        Add GB details which will be annotated using PLDO

        Parameters
        ----------
        gb_dict: dict
            dict containing details about the grain boundary

        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        #mark that the structure has a defect

        plane_defect_01 = BNode(name)
        self.add((self.material, CMSO.hasDefect, plane_defect_01))
        
        if gb_dict["GBType"] is None:
            self.add((plane_defect_01, RDF.type, PLDO.GrainBoundary))
        elif gb_dict["GBType"] == "Twist":
            self.add((plane_defect_01, RDF.type, PLDO.TwistBoundary))
        elif gb_dict["GBType"] == "Tilt":
            self.add((plane_defect_01, RDF.type, PLDO.TiltBoundary))
        elif gb_dict["GBType"] == "Symmetric Tilt":
            self.add((plane_defect_01, RDF.type, PLDO.SymmetricTiltBoundary))
        elif gb_dict["GBType"] == "Mixed":
            self.add((plane_defect_01, RDF.type, PLDO.MixedBoundary))
        self.add((plane_defect_01, PLDO.hasSigmaValue, Literal(gb_dict["sigma"], datatype=XSD.integer)))
        
        #now mark that the defect is GB
        #uname = None
        #if name is not None:
        #    uname = f'{name}GrainBoundaryPlane'
        #gb_plane_01 = BNode(uname)
        self.add((plane_defect_01, PLDO.hasGBPlane, Literal(gb_dict["GBPlane"], 
                                                             datatype=XSD.string)))
        #self.add((gb_plane_01, RDF.type, PLDO.GrainBoundaryPlane))
        #self.add((gb_plane_01, PLDO.hasMillerIndices, Literal(gb_dict["GBPlane"], 
        #                                                     datatype=XSD.string)))
        
        #uname = None
        #if name is not None:
        #    uname = f'{name}RotationAxis'
        #rotation_axis_01 = BNode(uname)
        self.add((plane_defect_01, PLDO.hasRotationAxis, Literal(gb_dict["RotationAxis"], 
                                                             datatype=XSD.string)))
        #self.add((rotation_axis_01, RDF.type, PLDO.RotationAxis))
        #self.add((rotation_axis_01, PLDO.hasComponentX, Literal(gb_dict["RotationAxis"][0], datatype=XSD.float)))
        #self.add((rotation_axis_01, PLDO.hasComponentY, Literal(gb_dict["RotationAxis"][1], datatype=XSD.float)))
        #self.add((rotation_axis_01, PLDO.hasComponentZ, Literal(gb_dict["RotationAxis"][2], datatype=XSD.float)))

        #uname = None
        #if name is not None:
        #    uname = f'{name}MisorientationAngle'
        #misorientation_angle_01 = BNode(uname)
        self.add((plane_defect_01, PLDO.hasMisorientationAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))
        #self.add((misorientation_angle_01, RDF.type, PLDO.MisorientationAngle))
        #self.add((misorientation_angle_01, PLDO.hasAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))    
    
    def add_vacancy(self, concentration, number=None, name=None):
        """
        Add Vacancy details which will be annotated by PODO

        Parameters
        ----------
        concentration: float
            vacancy concentration, value should be between 0-1

        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        vacancy_01 = BNode(name)
        self.add((self.material, CMSO.hasDefect, vacancy_01))
        self.add((vacancy_01, RDF.type, PODO.Vacancy))
        self.add((self.simulation_cell, PODO.hasVacancyConcentration, Literal(concentration, datatype=XSD.float)))
        if number is not None:
            self.add((self.simulation_cell, PODO.hasNumberOfVacancies, Literal(number, datatype=XSD.integer)))
        #if vacancy is added, atoms have to be deleted from the existing record!
        #this is indeed a tricky item


    def add_calculated_quantity(self, propertyname, value, unit=None, sample=None):
        prop = BNode()
        if sample is None:
            sample = self.sample
        self.add((sample, CMSO.hasCalculatedProperty, prop))
        self.add((prop, RDF.type, CMSO.CalculatedProperty))
        self.add((prop, RDFS.label, Literal(propertyname)))
        self.add((prop, CMSO.hasValue, Literal(value)))
        if unit is not None:
            self.add((prop, CMSO.hasUnit, URIRef(f'https://qudt.org/2.1/vocab/unit#{unit}')))


    def inspect_sample(self, sample=None):
        if sample is None:
            sample = self.sample
        natoms = self.graph.value(sample, CMSO.hasNumberOfAtoms).toPython()
        material = list([k[2] for k in self.graph.triples((sample, CMSO.hasMaterial, None))])[0]
        defects = list([k[2] for k in self.graph.triples((material, CMSO.hasDefect, None))])
        composition = list([k[2].toPython() for k in self.graph.triples((material, CMSO.hasElementRatio, None))])
        crystalstructure = self.graph.value(material, CMSO.hasStructure)
        spacegroupsymbol = self.graph.value(crystalstructure, CMSO.hasSpaceGroupSymbol).toPython()

        lattice = self.graph.value(sample, CMSO.hasNumberOfAtoms).toPython()
        defect_types = list([self.graph.value(d, RDF.type).toPython() for d in defects])
        prop_nodes = list([k[2] for k in self.graph.triples((sample, CMSO.hasCalculatedProperty, None))])
        props = list([self.graph.value(prop_node, RDFS.label) for prop_node in prop_nodes])
        propvals = list([self.graph.value(d, CMSO.hasValue).toPython() for d in prop_nodes])
        units = list([self.graph.value(d, CMSO.hasUnit).toPython() for d in prop_nodes])
        st = []
        st.append(f'Sample with {natoms} atoms.\n')
        st.append("Material:\n")
        st.append(" ".join(composition))
        st.append("\n")
        st.append(f'Space Group symbol: {spacegroupsymbol}\n')
        if len(defect_types) > 0:
            st.append('With defects:\n')
            for d in defect_types:
                st.append(f'{d}\n')
        if len(props) > 0:
            st.append('With calculated properties:\n')
            for x in range(len(props)):
                st.append(f'{props[x]} with value: {propvals[x]} and unit: {units[x]}\n')

        return " ".join(st)

    def visualize(self, *args, **kwargs):
        """
        Vosualise the RDF tree of the Graph

        Parameters
        ----------
        backend: string, {'ipycytoscape', 'graphviz'}
            Chooses the backend with which the graph will be plotted. ipycytoscape provides an interactive, 
            but slow visualisation, whereas graphviz provides a non-interactive fast visualisation.

        edge_color: string
            Edge color of the boxes

        styledict: dict
            If provided, allows customisation of color and other properties.

        graph_attr: dict
            further attributes that allow customisation of graphs

        layoutname: string
            name of the layout for graph

        Returns
        -------

        Notes
        -----
        styledict has the following options. Refer to graphviz and ipycytoscape
        documentation for more details
        BNode:  
          color:  
          shape:  
          style:  
        URIRef:  
          color:  
          shape:  
          style:  
        Literal:  
          color:  
          shape:  
          style:         
        """
        self.visualise(*args, **kwargs)
        
    def visualise(self,
                  backend='ipycytoscape',
                  edge_color="#37474F",
                  styledict=None, 
                  graph_attr ={'rankdir': 'BT'},
                  layoutname='cola'):
        """
        Vosualise the RDF tree of the Graph

        Parameters
        ----------
        backend: string, {'ipycytoscape', 'graphviz'}
            Chooses the backend with which the graph will be plotted. ipycytoscape provides an interactive, 
            but slow visualisation, whereas graphviz provides a non-interactive fast visualisation.

        edge_color: string
            Edge color of the boxes

        styledict: dict
            If provided, allows customisation of color and other properties.

        graph_attr: dict
            further attributes that allow customisation of graphs

        layoutname: string
            name of the layout for graph

        Returns
        -------

        Notes
        -----
        styledict has the following options. Refer to graphviz and ipycytoscape
        documentation for more details
        BNode:  
          color:  
          shape:  
          style:  
        URIRef:  
          color:  
          shape:  
          style:  
        Literal:  
          color:  
          shape:  
          style:         
        """
        
        sdict = defstyledict.copy()
        if styledict is not None:
            sdict = _replace_keys(sdict, styledict)
        return visualize_graph(self.graph, 
                               backend=backend,
                               edge_color=edge_color,
                               styledict=sdict, 
                               graph_attr=graph_attr,
                               layoutname=layoutname)
    
    
    def write(self, filename, format="json-ld"):
        """
        Write the serialised version of the graph to a file

        Parameters
        ----------
        filename: string
            name of output file

        format: string, {'turtle', 'xml', 'json-ld', 'ntriples', 'n3'}
            output format to be written to 

        Returns
        -------
        None
        """

        with open(filename, "w") as fout:
            fout.write(self.graph.serialize(format=format))
            
        
    def to_file(self, sample, filename=None, format="lammps-dump"):
        """
        Save a given sample to a file

        Parameters
        ----------
        sample
            ID of the sample

        filename: string
            name of output file

        format: string, {"lammps-dump","lammps-data", "poscar"}

        Returns
        -------
        None
        """

        if filename is None:
            filename = os.path.join(os.getcwd(), "out")
        sys = self.get_system_from_sample(sample)
        if format=="ase":
            return sys.to_ase()
        elif format=='poscar':
            asesys = sys.to_ase()
            write(filename, asesys, format="vasp")
        else:
            #asesys = sys.to_ase()
            #write(filename, asesys, format=format)
            sys.to_file(filename, format=format)
    
    
    def query(self, inquery):
        """
        Query the graph using SPARQL

        Parameters
        ----------
        inquery: string
            SPARQL query to be executed

        Returns
        -------
        res: pandas DataFrame
            pandas dataframe results
        """
        res = self.graph.query(inquery)
        if res is not None:
            return pd.DataFrame(res)
        raise ValueError("SPARQL query returned None")



    def query_sample(self, target_property, value, return_query=False):
        """
        Query the Graph for a sample that has the given `value` for the given `target_property`

        Parameters
        ----------
        target_property: string
            The target property can be any Ontology data property

        value: number, string, or list of either
            The value of the target property to be queried

        return_query: bool, optional
            If True, return the SPARQL query 

        Returns
        -------
        res: list
            list of queried samples

        query: string
            only returned if `return_query` is True

        Notes
        -----
        """

        query = self._query_graph.formulate_query(target_property, value)
        res = self.graph.query(query)
        res = [r for r in res]
        if return_query:
            return res, query
        return res
    
    #################################
    # Methods to interact with sample
    #################################
    @property
    def n_samples(self):
        """
        Number of samples in the Graph
        """

        return len([x for x in self.graph.triples((None, RDF.type, CMSO.AtomicScaleSample))])
    
    @property
    def samples(self):
        """
        Returns a list of all Samples in the graph
        """

        return [x[0] for x in self.graph.triples((None, RDF.type, CMSO.AtomicScaleSample))]
        
    def iterate_graph(self, item, create_new_graph=False):
        if create_new_graph:
            self.sgraph = RDFGraph()
        triples = list(self.graph.triples((item, None, None)))
        for triple in triples:
            self.sgraph.graph.add(triple)
            self.iterate_graph(triple[2])
    
    def get_sample(self, sample, no_atoms=False):
        """
        Get the Sample as an RDFGraph

        Parameters
        ----------
        sample: string
            sample id

        no_atoms: bool, optional
            if True, returns the number of atoms in the sample

        Returns
        -------
        sgraph: :py:class:`RDFGraph`
            the RDFGraph of the queried sample

        """

        self.iterate_graph(sample, create_new_graph=True)
        if no_atoms:
            na = self.sgraph.graph.value(sample, CMSO.hasNumberOfAtoms).toPython()
            return self.sgraph, na
        return self.sgraph
        
    def get_system_from_sample(self, sample):
        """
        Get a pyscal :py:class:`pyscal.core.System` from the selected sample

        Parameters
        ----------
        sample: string
            sample id

        Returns
        -------
        system: :py:class:`pyscal.core.System`
            corresponding system
        """

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
