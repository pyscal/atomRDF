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
import yaml
import uuid
import json
import shutil
import tarfile
import pyscal_rdf.json_io as json_io

from pyscal_rdf.visualize import visualize_graph
from pyscal_rdf.network.network import OntologyNetwork
from pyscal_rdf.network.ontology import read_ontology
from pyscal_rdf.rdfsystem import System
import pyscal_rdf.properties as prp
#from pyscal3.core import System
from pyscal3.atoms import Atoms

CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("http://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")

#read element data file
file_location = os.path.dirname(__file__).split('/')
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__),  'data/element.yml')
with open(file_location, 'r') as fin:
    element_indetifiers = yaml.safe_load(fin)


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

def _setup_structure_store(structure_store):
    if structure_store is None:
        structure_store = os.path.join(os.getcwd(), 'rdf_structure_store')
    if not os.path.exists(structure_store):
        os.mkdir(structure_store)
    return structure_store

class RDFGraph:
    def __init__(self, graph_file=None, 
        store="Memory", 
        store_file=None,
        identifier="http://default_graph",
        ontology=None,
        structure_store=None):
        
        self.store_file = store_file
        self.structure_store = structure_store


        if store == "Memory":
            self.graph = Graph(store="Memory", identifier=identifier)

                
        elif store=="SQLAlchemy":
            if store_file is None:
                raise ValueError("store file is needed if store is not memory")
            self.graph = Graph(store="SQLAlchemy", identifier=identifier)
            uri = Literal(f"sqlite:///{store_file}")
            self.graph.open(uri, create=True)        
        else:
            raise ValueError("store should be pyiron_project, SQLAlchemy, or Memory")
        
        #start the storage
        self.structure_store = _setup_structure_store(self.structure_store)

        #start binding
        self.graph.bind("cmso", CMSO)
        self.graph.bind("pldo", PLDO)
        
        if graph_file is not None:
            if os.path.exists(graph_file):
                self.graph.parse(graph_file)
        
        self.sample = None
        self.material = None
        self.sysdict = None
        self.sgraph = None
        if ontology is None:
            ontology = read_ontology()
        self.ontology = ontology
        self.terms = self.ontology.terms
        self._atom_ids = None
        self.store = store

    
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
            self.system = structure
        elif os.path.exists(structure):
            sys = System(structure, format=format)
            #self.sysdict = convert_to_dict(sys)
            self.system = sys
        
    def add(self, triple):
        if str(triple[2].toPython()) != 'None':
            self.graph.add(triple)
        
    def add_structure_to_graph(self, 
            structure, 
            names=False, 
            name_index=None, 
            format=None):
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

        Notes
        -----
        BNodes, or relational nodes will be avoided as much as possible so that merging of datasets would be possible.
        Instead URIref containers will be made use of. This makes the `names` and `name_index` parameters crucial.
        `names` parameter means that legible names starting with the string `Sample_x` would be used. `x` would ensure
        that there is conflict with the current database. However, they do not ensure there is no conflicts when various
        graphs are merged together. Hence this value is recommended only for simple, demonstration cases.

        If `names` are False, unique ids are generated which would be id of the sample. These ids use the python `uuid` module
        and therefore ensures that the names are always unique.
        """
        
        self.process_structure(structure, format=format)
        
        if names:
            if name_index is None:
                name_index = self.n_samples + 1
            self._name = f'sample:{name_index}'
        else:
            self._name = f'sample:{str(uuid.uuid4())}'

        self.create_graph()
        structure.sample = self.sample
        #structure._atom_ids = copy.copy(self._atom_ids)
        structure.graph = self
    
    def create_graph(self):
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
        self.add_sample()
        self.add_material()
        self.add_chemical_composition()
        self.add_simulation_cell()
        self.add_simulation_cell_properties()
        self.add_crystal_structure()
        self.add_space_group()
        self.add_unit_cell()
        self.add_lattice_properties()
        self.add_atoms()

        #extra triples
        self.add((CMSO.SimulationCellLength, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.LatticeParameter, RDFS.subClassOf, CMSO.Length))
        self.add((CMSO.Length, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM")))
        
        self.add((CMSO.SimulationCellAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.LatticeAngle, RDFS.subClassOf, CMSO.Angle))
        self.add((CMSO.Angle, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/DEG")))
        
        self.add((CMSO.LatticeVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.SimulationCellVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.PositionVector, RDFS.subClassOf, CMSO.Vector))
        self.add((CMSO.Vector, CMSO.hasUnit, URIRef("http://qudt.org/vocab/unit/ANGSTROM")))
        
        
    def add_sample(self):
        """
        Add a CMSO Sample object

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        sample = URIRef(f'{self._name}')
        self.add((sample, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample
    
    def add_material(self):
        """
        Add a CMSO Material object

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        material = URIRef(f'{self._name}_Material')
        self.add((self.sample, CMSO.hasMaterial, material))
        self.add((material, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material
    
    def add_chemical_composition(self):
        """
        Add chemical composition

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        composition = self.system.schema.material.element_ratio()

        chemical_species = URIRef(f'{self._name}_ChemicalSpecies')
        self.add((self.sample, CMSO.hasSpecies, chemical_species))
        self.add((chemical_species, RDF.type, CMSO.ChemicalSpecies))

        for e, r in composition.items():
            if e in element_indetifiers.keys():
                element = URIRef(element_indetifiers[e])
                self.add((chemical_species, CMSO.hasElement, element))
                self.add((element, RDF.type, CMSO.Element))
                self.add((element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string)))
                self.add((element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float)))
    
    def add_simulation_cell(self):
        """
        Add a CMSO SimulationCell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        simulation_cell = URIRef(f'{self._name}_SimulationCell')
        self.add((self.sample, CMSO.hasSimulationCell, simulation_cell))
        self.add((simulation_cell, RDF.type, CMSO.SimulationCell))
        self.add((simulation_cell, CMSO.hasVolume, 
            Literal(np.round(self.system.schema.simulation_cell.volume(), decimals=2), 
                datatype=XSD.float)))
        self.add((self.sample, CMSO.hasNumberOfAtoms, 
            Literal(self.system.schema.simulation_cell.number_of_atoms(), 
                datatype=XSD.integer)))
        self.simulation_cell = simulation_cell
        
    
    def add_simulation_cell_properties(self):
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
        simulation_cell_length = URIRef(f'{self._name}_SimulationCellLength')
        self.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length))
        data = self.system.schema.simulation_cell.length()
        self.add((simulation_cell_length, RDF.type, CMSO.SimulationCellLength))
        self.add((simulation_cell_length, CMSO.hasLength_x, Literal(data[0], datatype=XSD.float)))
        self.add((simulation_cell_length, CMSO.hasLength_y, Literal(data[1], datatype=XSD.float)))
        self.add((simulation_cell_length, CMSO.hasLength_z, Literal(data[2], datatype=XSD.float)))
        
        simulation_cell_vector_01 = URIRef(f'{self._name}_SimulationCellVector_1')
        data = self.system.schema.simulation_cell.vector()
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        self.add((simulation_cell_vector_01, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_x, Literal(data[0][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_y, Literal(data[0][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_01, CMSO.hasComponent_z, Literal(data[0][2], datatype=XSD.float)))
        
        simulation_cell_vector_02 = URIRef(f'{self._name}_SimulationCellVector_2')
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        self.add((simulation_cell_vector_02, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_x, Literal(data[1][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_y, Literal(data[1][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_02, CMSO.hasComponent_z, Literal(data[1][2], datatype=XSD.float)))
        
        simulation_cell_vector_03 = URIRef(f'{self._name}_SimulationCellVector_3')
        self.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        self.add((simulation_cell_vector_03, RDF.type, CMSO.SimulationCellVector))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_x, Literal(data[2][0], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_y, Literal(data[2][1], datatype=XSD.float)))
        self.add((simulation_cell_vector_03, CMSO.hasComponent_z, Literal(data[2][2], datatype=XSD.float)))
        
        simulation_cell_angle = URIRef(f'{self._name}_SimulationCellAngle')
        data = self.system.schema.simulation_cell.angle()
        self.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle))
        self.add((simulation_cell_angle, RDF.type, CMSO.SimulationCellAngle))
        self.add((simulation_cell_angle, CMSO.hasAngle_alpha, Literal(data[0], datatype=XSD.float)))
        self.add((simulation_cell_angle, CMSO.hasAngle_beta, Literal(data[1], datatype=XSD.float)))
        self.add((simulation_cell_angle, CMSO.hasAngle_gamma, Literal(data[2], datatype=XSD.float)))
        
    
    def add_crystal_structure(self):
        """
        Add a CMSO Crystal Structure

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        crystal_structure = URIRef(f'{self._name}_CrystalStructure')
        self.add((self.material, CMSO.hasStructure, crystal_structure))
        self.add((crystal_structure, RDF.type, CMSO.CrystalStructure))    
        self.add((crystal_structure, CMSO.hasAltName, 
            Literal(self.system.schema.material.crystal_structure.name(), 
                datatype=XSD.string)))
        self.crystal_structure = crystal_structure
        
    def add_space_group(self):
        """
        Add a CMSO Space Group

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        space_group = URIRef(f'{self._name}_SpaceGroup')
        self.add((self.crystal_structure, CMSO.hasSpaceGroup, space_group))
        self.add((space_group, CMSO.hasSpaceGroupSymbol, 
            Literal(self.system.schema.material.crystal_structure.spacegroup_symbol(), 
                datatype=XSD.string)))
        self.add((space_group, CMSO.hasSpaceGroupNumber, 
            Literal(self.system.schema.material.crystal_structure.spacegroup_number(), 
                datatype=XSD.integer)))
    
            
    def add_unit_cell(self):
        """
        Add a CMSO Unit Cell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        unit_cell = URIRef(f'{self._name}_UnitCell')
        self.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell))
        self.add((unit_cell, RDF.type, CMSO.UnitCell))
        self.unit_cell = unit_cell
        
        #add bravais lattice
        bv = self.system.schema.material.crystal_structure.unit_cell.bravais_lattice()
        if bv is not None:
            bv = URIRef(bv)
            self.add((self.unit_cell, CMSO.hasBravaisLattice, bv))
        
    def add_lattice_properties(self):
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
        data = self.system.schema.material.crystal_structure.unit_cell.lattice_parameter()
        lattice_parameter = URIRef(f'{self._name}_LatticeParameter')
        self.add((self.unit_cell, CMSO.hasLatticeParamter, lattice_parameter))
        self.add((lattice_parameter, RDF.type, CMSO.LatticeParameter))
        self.add((lattice_parameter, CMSO.hasLength_x, Literal(data[0], datatype=XSD.float)))
        self.add((lattice_parameter, CMSO.hasLength_y, Literal(data[1], datatype=XSD.float)))
        self.add((lattice_parameter, CMSO.hasLength_z, Literal(data[2], datatype=XSD.float)))
        
        lattice_angle = URIRef(f'{self._name}_LatticeAngle')
        data = self.system.schema.material.crystal_structure.unit_cell.angle()
        self.add((self.unit_cell, CMSO.hasAngle, lattice_angle))
        self.add((lattice_angle, RDF.type, CMSO.LatticeAngle))
        self.add((lattice_angle, CMSO.hasAngle_alpha, Literal(data[0], datatype=XSD.float)))
        self.add((lattice_angle, CMSO.hasAngle_beta, Literal(data[1], datatype=XSD.float)))
        self.add((lattice_angle, CMSO.hasAngle_gamma, Literal(data[2], datatype=XSD.float)))        


    def _save_atom_attributes(self, position_identifier, species_identifier):
        #if self.store == 'pyiron':
        #    pass
        #else:
        #    #this is the file based store system
        datadict = {
            position_identifier:{
                "value": self.system.schema.atom_attribute.position(),
                "label": "position", 
            },
            species_identifier:{
                "value": self.system.schema.atom_attribute.species(),
                "label": "species", 
            },
        }
        outfile = os.path.join(self.structure_store, str(self._name).split(':')[-1])
        json_io.write_file(outfile,  datadict)
        return os.path.relpath(outfile+'.json')


    def add_atoms(self):
        """
        Add Atoms including their species and positions

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------

        Notes
        -----
        Note that for the moment, we will dump the structures in a given folder,
        maybe this could be input from the Job class directly
        """
        #now we write out file
        position_identifier = str(uuid.uuid4())
        species_identifier = str(uuid.uuid4())

        outfile = self._save_atom_attributes(position_identifier, species_identifier)

        if "positions" in self.system.atoms.keys():
            position = URIRef(f'{self._name}_Position')
            self.add((self.sample, CMSO.hasAttribute, position))
            self.add((position, RDF.type, CMSO.AtomAttribute))
            self.add((position, CMSO.hasName, Literal('Position', datatype=XSD.string)))
            self.add((position, CMSO.hasIdentifier, Literal(position_identifier, datatype=XSD.string)))            
            self.add((position, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))

        if "species" in self.system.atoms.keys():
            species = URIRef(f'{self._name}_Species')
            self.add((self.sample, CMSO.hasAttribute, species))
            self.add((species, RDF.type, CMSO.AtomAttribute))
            self.add((species, CMSO.hasName, Literal('Species', datatype=XSD.string)))
            self.add((species, CMSO.hasIdentifier, Literal(species_identifier, datatype=XSD.string)))            
            self.add((species, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))

        #if "velocities" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Velocity'
        #    velocity = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, velocity))
        #    self.add((velocity, RDF.type, CMSO.AtomAttribute))
        #    self.add((velocity, CMSO.hasName, Literal('Velocity', data_type=XSD.string)))
        #    velocity_identifier = uuid.uuid4()
        #    self.add((velocity, CMSO.hasIdentifier, Literal(velocity_identifier, datatype=XSD.string)))            

        #if "forces" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Force'  
        #    force = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, force))
        #    self.add((force, RDF.type, CMSO.AtomAttribute))
        #    self.add((force, CMSO.hasName, Literal('Force', data_type=XSD.string)))
        #    force_identifier = uuid.uuid4()
        #    self.add((force, CMSO.hasIdentifier, Literal(force_identifier, datatype=XSD.string)))            


    
    
    def add_gb(self, gb_dict):
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
        if gb_dict["GBType"] is None:
            plane_defect = URIRef(f'{self._name}_GrainBoundary')
            self.add((plane_defect, RDF.type, PLDO.GrainBoundary))
        
        elif gb_dict["GBType"] == "Twist":
            plane_defect = URIRef(f'{self._name}_TwistGrainBoundary')
            self.add((plane_defect, RDF.type, PLDO.TwistGrainBoundary))
        
        elif gb_dict["GBType"] == "Tilt":
            plane_defect = URIRef(f'{self._name}_TiltGrainBoundary')
            self.add((plane_defect, RDF.type, PLDO.TiltGrainBoundary))
        
        elif gb_dict["GBType"] == "Symmetric Tilt":
            plane_defect = URIRef(f'{self._name}_SymmetricalTiltGrainBoundary')
            self.add((plane_defect, RDF.type, PLDO.SymmetricalTiltGrainBoundary))
        
        elif gb_dict["GBType"] == "Mixed":
            plane_defect = URIRef(f'{self._name}_MixedGrainBoundary')
            self.add((plane_defect, RDF.type, PLDO.MixedGrainBoundary))
        
        self.add((self.material, CMSO.hasDefect, plane_defect))
        self.add((plane_defect, PLDO.hasSigmaValue, Literal(gb_dict["sigma"], datatype=XSD.integer)))
        
        #now mark that the defect is GB
        #uname = None
        #if name is not None:
        #    uname = f'{name}GrainBoundaryPlane'
        #gb_plane_01 = BNode(uname)
        self.add((plane_defect, PLDO.hasGBPlane, Literal(gb_dict["GBPlane"], 
                                                             datatype=XSD.string)))
        #self.add((gb_plane_01, RDF.type, PLDO.GrainBoundaryPlane))
        #self.add((gb_plane_01, PLDO.hasMillerIndices, Literal(gb_dict["GBPlane"], 
        #                                                     datatype=XSD.string)))
        
        #uname = None
        #if name is not None:
        #    uname = f'{name}RotationAxis'
        #rotation_axis_01 = BNode(uname)
        self.add((plane_defect, PLDO.hasRotationAxis, Literal(gb_dict["RotationAxis"], 
                                                             datatype=XSD.string)))
        #self.add((rotation_axis_01, RDF.type, PLDO.RotationAxis))
        #self.add((rotation_axis_01, PLDO.hasComponentX, Literal(gb_dict["RotationAxis"][0], datatype=XSD.float)))
        #self.add((rotation_axis_01, PLDO.hasComponentY, Literal(gb_dict["RotationAxis"][1], datatype=XSD.float)))
        #self.add((rotation_axis_01, PLDO.hasComponentZ, Literal(gb_dict["RotationAxis"][2], datatype=XSD.float)))

        #uname = None
        #if name is not None:
        #    uname = f'{name}MisorientationAngle'
        #misorientation_angle_01 = BNode(uname)
        self.add((plane_defect, PLDO.hasMisorientationAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))
        #self.add((misorientation_angle_01, RDF.type, PLDO.MisorientationAngle))
        #self.add((misorientation_angle_01, PLDO.hasAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))    
    
    def add_vacancy(self, concentration, number=None):
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

        vacancy = URIRef(f'{self._name}_Vacancy')
        self.add((self.material, CMSO.hasDefect, vacancy))
        self.add((vacancy, RDF.type, PODO.Vacancy))
        self.add((self.simulation_cell, PODO.hasVacancyConcentration, Literal(concentration, datatype=XSD.float)))
        if number is not None:
            self.add((self.simulation_cell, PODO.hasNumberOfVacancies, Literal(number, datatype=XSD.integer)))
        #if vacancy is added, atoms have to be deleted from the existing record!
        #this is indeed a tricky item


    def add_calculated_quantity(self, propertyname, value, unit=None, sample=None):
        prop = URIRef(f'{self._name}_{propertyname}')
        if sample is None:
            sample = self.sample
        self.add((sample, CMSO.hasCalculatedProperty, prop))
        self.add((prop, RDF.type, CMSO.CalculatedProperty))
        self.add((prop, RDFS.label, Literal(propertyname)))
        self.add((prop, CMSO.hasValue, Literal(value)))
        if unit is not None:
            self.add((prop, CMSO.hasUnit, URIRef(f'http://qudt.org/vocab/unit/{unit}')))


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
                  layoutname='cola',
                  hide_types=False,
                  workflow_view=False):
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
                               layoutname=layoutname,
                               hide_types=hide_types,
                               workflow_view=workflow_view)
    
    
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
    
    def archive(self, package_name, format='turtle', compress=True):
        """
        Publish a dataset from graph including per atom quantities
        """
        #first step make a folder
        if os.path.exists(package_name):
            raise ValueError(f'{package_name} already exists')
        if compress:
            if os.path.exists(f'{package_name}.tar.gz'):
                raise ValueError(f'{package_name} tarball already exists')
        
        os.mkdir(package_name)
        structure_store = f'{package_name}/rdf_structure_store' 
        os.mkdir(structure_store)

        #now go through each sample, and copy the file, at the same time fix the paths
        for sample in self.samples:
            filepath = self.graph.value(URIRef(f'{sample}_Position'), CMSO.hasPath).toPython()
            shutil.copy(filepath, structure_store)
            
            #now we have to remove the old path, and fix new
            for val in ['Position', 'Species']:
                self.graph.remove((URIRef(f'{sample}_{val}'), CMSO.hasPath, None))
            
                #assign corrected path
                new_relpath = "/".join(['rdf_structure_store', filepath.split('/')[-1]])
                self.graph.add((URIRef(f'{sample}_{val}'), CMSO.hasPath, Literal(new_relpath, datatype=XSD.string)))

        triple_file = os.path.join(package_name, 'triples')
        self.write(triple_file, format=format)

        if compress:
            with tarfile.open(f'{package_name}.tar.gz', "w:gz") as tar:
                tar.add(package_name, arcname=os.path.basename(package_name))
            shutil.rmtree(package_name)


    @classmethod
    def unarchive(cls, package_name, compress=True, 
        store="Memory", 
        store_file=None,
        identifier="http://default_graph",
        ontology=None):
        if compress:
            package_base_name = ".".join(package_name.split(".")[:-2])
            with tarfile.open(package_name) as fin: 
                fin.extractall(".")
            #os.remove(package_name)
            #copy things out

        return cls(store=store, store_file=store_file,
            identifier=identifier, 
            graph_file=f'{package_base_name}/triples', 
            structure_store=f'{package_base_name}/rdf_structure_store',
            ontology=ontology)
    
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
            for line in inquery.split('\n'):
                if 'SELECT DISTINCT' in line:
                    break
            labels = [x[1:] for x in line.split()[2:]]
            return pd.DataFrame(res, columns=labels)
        raise ValueError("SPARQL query returned None")

    def auto_query(self, source, destination, 
        condition=None, 
        return_query=False, 
        enforce_types=None):

        if enforce_types is None:
            for val in [True, False]:
                query = self.ontology.create_query(source, destination, 
                    condition=condition, enforce_types=val)
                if return_query:
                    return query
                res = self.query(query)
                if len(res) != 0:
                    return res
        else:
            query = self.ontology.create_query(source, destination, 
                condition=condition, enforce_types=val)
            if return_query:
                return query
            res = self.query(query)

        return res    

    #################################
    # Methods to interact with sample
    #################################
    def query_sample(self, destination, condition=None, return_query=False, enforce_types=True):
        return self.auto_query(self.ontology.terms.cmso.AtomicScaleSample, destination,
            condition=condition, return_query=return_query, enforce_types=enforce_types)

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
        filepath = self.graph.value(URIRef(f'{sample}_Position'), CMSO.hasPath).toPython()
        position_identifier = self.graph.value(URIRef(f'{sample}_Position'), CMSO.hasIdentifier).toPython()
        species_identifier = self.graph.value(URIRef(f'{sample}_Species'), CMSO.hasIdentifier).toPython()

        #open the file for reading
        with open(filepath, 'r') as fin:
            data = json.load(fin)
            positions = data[position_identifier]['value']
            species = data[species_identifier]['value']

        atoms = {"positions": positions, "species": species}
        at = Atoms()
        at.from_dict(atoms)
        sys = System()
        sys.box = cell_vectors
        sys.atoms = at       
        return sys

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
            return sys.write.ase()
        elif format=='poscar':
            asesys = sys.write.ase()
            write(filename, asesys, format="vasp")
        else:
            asesys = sys.write.ase()
            write(filename, asesys, format=format)
