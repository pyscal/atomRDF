"""
StructureGraph is the central object in pyscal_rdf which combines all the functionality
of :py:class:`pyscal_rdf.graph.RDFGraph` along with easy structural creation routines.
"""
import numpy as np
import copy
from functools import partial, update_wrapper
import os
import yaml
import uuid
import json
import shutil
import tarfile

import pyscal_rdf.json_io as json_io

import pyscal3.structure_creator as pcs
from pyscal3.grain_boundary import GrainBoundary
from pyscal3.atoms import AttrSetter
import pyscal_rdf.properties as prp
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict

from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("http://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")

#read element data file
file_location = os.path.dirname(__file__).split('/')
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__),  'data/element.yml')
with open(file_location, 'r') as fin:
    element_indetifiers = yaml.safe_load(fin)

def _make_crystal(structure, 
    lattice_constant = 1.00, 
    repetitions = None, 
    ca_ratio = 1.633, 
    noise = 0, 
    element=None,
    primitive=False,
    graph=None,
    names=True):
    
    atoms, box, sdict = pcs.make_crystal(structure, 
        lattice_constant=lattice_constant,
        repetitions=repetitions, 
        ca_ratio=ca_ratio,
        noise=noise, 
        element=element, 
        return_structure_dict=True,
        primitive=primitive)
    
    s = System(graph=graph, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    s.to_graph()
    return s

def _make_general_lattice(positions,
    types, 
    box,
    lattice_constant = 1.00, 
    repetitions = None, 
    noise = 0,
    element=None,
    graph=None,
    names=True):

    atoms, box, sdict = pcs.general_lattice(positions,
        types,
        box,
        lattice_constant=lattice_constant,
        repetitions=repetitions,
        noise=noise,
        element=element,
        return_structure_dict=True)
    s = System(graph=graph, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = 'custom'
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    s.to_graph()

    return s

def _make_grain_boundary(axis, 
    sigma, gb_plane,
    structure = None,
    element = None, 
    lattice_constant = 1,
    repetitions = (1,1,1),
    overlap=0.0,
    graph=None,
    names=True):

    gb = GrainBoundary()
    gb.create_grain_boundary(axis=axis, sigma=sigma, 
                             gb_plane=gb_plane)

    if structure is not None:
        atoms, box, sdict = gb.populate_grain_boundary(structure, 
                                        repetitions = repetitions,
                                        lattice_parameter = lattice_constant,
                                        overlap=overlap)
    elif element is not None:
        atoms, box, sdict = gb.populate_grain_boundary(element, 
                                        repetitions=repetitions,
                                        overlap=overlap)
    s = System(graph=graph, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    s.to_graph()
    gb_dict = {"GBPlane": " ".join(np.array(gb_plane).astype(str)),
              "RotationAxis": axis,
              "MisorientationAngle": gb.theta,
              "GBType": gb.find_gb_character(),
              "sigma": gb.sigma,
              }
    s.add_gb(gb_dict)
    return s

def _read_structure(filename, 
                  format="lammps-dump",
                  graph=None, 
                  names=False,
                  species=None,
                  lattice=None,
                  lattice_constant=None,
                  basis_box=None,
                  basis_positions=None,
                  ):
    datadict = {}
    if lattice is not None:
        if lattice in structure_dict.keys():
            datadict = structure_dict[lattice]['conventional']
        datadict['lattice'] = lattice
    if lattice_constant is not None:
        datadict['lattice_constant'] = lattice_constant
    if basis_box is not None:
        datadict['box'] = basis_box
    if basis_positions is not None:
        datadict['positions'] = basis_positions

    s = System(filename, format=format, species=species, 
        graph=graph, names=names)
    s.lattice_properties = datadict
    s.to_graph()
    return s   



class System(pc.System):

    create = AttrSetter()
    #create.head = pcs
    mapdict = {}
    mapdict["lattice"] = {}
    for key in structure_dict.keys():
        mapdict["lattice"][key] = update_wrapper(partial(_make_crystal, key), 
            _make_crystal)
    mapdict["lattice"]["custom"] = _make_general_lattice

    mapdict["element"] = {}
    for key in element_dict.keys():
        mapdict["element"][key] = update_wrapper(partial(_make_crystal,
            element_dict[key]['structure'],
            lattice_constant=element_dict[key]['lattice_constant'],
            element = key), pcs.make_crystal)

    mapdict["defect"] = {}
    mapdict["defect"]["grain_boundary"] = _make_grain_boundary
    create._add_attribute(mapdict)

    read = AttrSetter()
    mapdict = {}
    mapdict['file'] = _read_structure
    mapdict['ase'] = update_wrapper(partial(_read_structure, format='ase'), _read_structure)
    read._add_attribute(mapdict)

    def __init__(self, filename = None, 
            format = "lammps-dump", 
            compressed = False, 
            customkeys = None,
            species = None,
            source=None,
            graph=None,
            names=True):
        
        super().__init__(filename = filename, 
            format = format, 
            compressed = compressed, 
            customkeys = customkeys,
            species = species)
        
        #this is the sample which will be stored
        self.sample = None
        #the graph object should also be attached
        #for post-processing of structures
        self.graph = graph
        self.names = names
        self._atom_ids = None
        if source is not None:
            self.__dict__.update(source.__dict__)

        #assign attributes
        self.schema = AttrSetter()
        mapdict = {
        "material": {
            "element_ratio": partial(prp.get_chemical_composition, self),
            "crystal_structure": {
                "name": partial(prp.get_crystal_structure_name, self),
                "spacegroup_symbol": partial(prp.get_spacegroup_symbol, self),
                "spacegroup_number": partial(prp.get_spacegroup_number, self),
                "unit_cell": {
                    "bravais_lattice": partial(prp.get_bravais_lattice, self),
                    "lattice_parameter": partial(prp.get_lattice_parameter, self),
                    "angle": partial(prp.get_lattice_angle, self),
                    },            
                },        
            },
        "simulation_cell": {
            "volume": partial(prp.get_cell_volume, self),
            "number_of_atoms": partial(prp.get_number_of_atoms, self),
            "length": partial(prp.get_simulation_cell_length, self),
            "vector": partial(prp.get_simulation_cell_vector, self),
            "angle": partial(prp.get_simulation_cell_angle, self),
            },
        "atom_attribute": {
            "position": partial(prp.get_position, self),
            "species": partial(prp.get_species, self),
            },
        }

        self.schema._add_attribute(mapdict)

    def delete(self, ids=None, indices=None, condition=None, selection=False):
        
        masks = self.atoms._generate_bool_list(ids=ids, indices=indices, condition=condition, selection=selection)
        delete_list = [masks[self.atoms["head"][x]] for x in range(self.atoms.ntotal)]
        delete_ids = [x for x in range(self.atoms.ntotal) if delete_list[x]]
        self.atoms._delete_atoms(delete_ids)

        if self.graph is not None:
            #first annotate graph
            val = len([x for x in masks if x])        
            c = (val/self.natoms)
            self.add_vacancy(c, number=val)
            #now we need to re-add atoms, so at to remove
            self.graph.graph.remove((self.sample, CMSO.hasNumberOfAtoms, None))
            self.graph.graph.add((self.sample, CMSO.hasNumberOfAtoms, Literal(self.natoms-val, datatype=XSD.integer)))
            #revamp composition
            #remove existing chem composution
            chemical_species = self.graph.graph.value(self.sample, CMSO.hasSpecies)
            #start by cleanly removing elements
            for s in self.graph.graph.triples((chemical_species, CMSO.hasElement, None)):
                element = s[2]
                self.graph.graph.remove((element, None, None))
            self.graph.graph.remove((chemical_species, None, None))
            self.graph.graph.remove((self.sample, CMSO.hasSpecies, None))
            
            #now recalculate and add it again
            composition = self.schema.material.element_ratio()

            chemical_species = URIRef(f'{self._name}_ChemicalSpecies')
            self.graph.graph.add((self.sample, CMSO.hasSpecies, chemical_species))
            self.graph.graph.add((chemical_species, RDF.type, CMSO.ChemicalSpecies))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = URIRef(element_indetifiers[e])
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add((element, RDF.type, CMSO.Element))
                    self.graph.add((element, CMSO.hasSymbol, Literal(e, datatype=XSD.string)))
                    self.graph.add((element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float)))

            #we also have to read in file and clean it up
            filepath = self.graph.graph.value(URIRef(f'{self.sample}_Position'), CMSO.hasPath).toPython()
            position_identifier = self.graph.graph.value(URIRef(f'{self.sample}_Position'), CMSO.hasIdentifier).toPython()
            species_identifier = self.graph.graph.value(URIRef(f'{self.sample}_Species'), CMSO.hasIdentifier).toPython()

            #open the file for reading
            with open(filepath, 'r') as fin:
                data = json.load(fin)
                positions = data[position_identifier]['value']
                species = data[species_identifier]['value']

            #clean up items
            positions = [pos for count, pos in enumerate(positions) if count not in delete_ids]
            species = [pos for count, pos in enumerate(species) if count not in delete_ids]

            datadict = {
                position_identifier:{
                    "value": positions,
                    "label": "position", 
                },
                species_identifier:{
                    "value": species,
                    "label": "species", 
                },
            }
            outfile = os.path.join(self.graph.structure_store, str(self._name).split(':')[-1])
            json_io.write_file(outfile,  datadict)


    def __delitem__(self, val):
        if isinstance(val, int):
            val = [val]
        #now the graph has to be updated accordingly
        self.delete(indices=list(val))

    def to_graph(self):
        if self.graph is None:
            return

        self._generate_name()
        self._add_sample()
        self._add_material()
        self._add_chemical_composition()
        self._add_simulation_cell()
        self._add_simulation_cell_properties()
        self._add_crystal_structure()
        self._add_space_group()
        self._add_unit_cell()
        self._add_lattice_properties()
        self._add_atoms()



    def _generate_name(self, name_index=None):
        if self.names:
            if name_index is None:
                name_index = self.graph.n_samples + 1
            self._name = f'sample:{name_index}'
        else:
            self._name = f'sample:{str(uuid.uuid4())}'

    def _add_sample(self):
        sample = URIRef(f'{self._name}')
        self.graph.add((sample, RDF.type, CMSO.AtomicScaleSample))
        self.sample = sample

    def _add_material(self):
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
        self.graph.add((self.sample, CMSO.hasMaterial, material))
        self.graph.add((material, RDF.type, CMSO.CrystallineMaterial))        
        self.material = material

    def _add_chemical_composition(self):
        """
        Add chemical composition

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        composition = self.schema.material.element_ratio()

        chemical_species = URIRef(f'{self._name}_ChemicalSpecies')
        self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))
        self.graph.add((chemical_species, RDF.type, CMSO.ChemicalSpecies))

        for e, r in composition.items():
            if e in element_indetifiers.keys():
                element = URIRef(element_indetifiers[e])
                self.graph.add((chemical_species, CMSO.hasElement, element))
                self.graph.add((element, RDF.type, CMSO.Element))
                self.graph.add((element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string)))
                self.graph.add((element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float)))

    def _add_simulation_cell(self):
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
        self.graph.add((self.sample, CMSO.hasSimulationCell, simulation_cell))
        self.graph.add((simulation_cell, RDF.type, CMSO.SimulationCell))
        self.graph.add((simulation_cell, CMSO.hasVolume, 
            Literal(np.round(self.schema.simulation_cell.volume(), decimals=2), 
                datatype=XSD.float)))
        self.graph.add((self.sample, CMSO.hasNumberOfAtoms, 
            Literal(self.schema.simulation_cell.number_of_atoms(), 
                datatype=XSD.integer)))
        self.simulation_cell = simulation_cell
        
    
    def _add_simulation_cell_properties(self):
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
        self.graph.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length))
        data = self.schema.simulation_cell.length()
        self.graph.add((simulation_cell_length, RDF.type, CMSO.SimulationCellLength))
        self.graph.add((simulation_cell_length, CMSO.hasLength_x, Literal(data[0], datatype=XSD.float)))
        self.graph.add((simulation_cell_length, CMSO.hasLength_y, Literal(data[1], datatype=XSD.float)))
        self.graph.add((simulation_cell_length, CMSO.hasLength_z, Literal(data[2], datatype=XSD.float)))
        
        simulation_cell_vector_01 = URIRef(f'{self._name}_SimulationCellVector_1')
        data = self.schema.simulation_cell.vector()
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        self.graph.add((simulation_cell_vector_01, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_x, Literal(data[0][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_y, Literal(data[0][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_01, CMSO.hasComponent_z, Literal(data[0][2], datatype=XSD.float)))
        
        simulation_cell_vector_02 = URIRef(f'{self._name}_SimulationCellVector_2')
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        self.graph.add((simulation_cell_vector_02, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_x, Literal(data[1][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_y, Literal(data[1][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_02, CMSO.hasComponent_z, Literal(data[1][2], datatype=XSD.float)))
        
        simulation_cell_vector_03 = URIRef(f'{self._name}_SimulationCellVector_3')
        self.graph.add((self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        self.graph.add((simulation_cell_vector_03, RDF.type, CMSO.SimulationCellVector))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_x, Literal(data[2][0], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_y, Literal(data[2][1], datatype=XSD.float)))
        self.graph.add((simulation_cell_vector_03, CMSO.hasComponent_z, Literal(data[2][2], datatype=XSD.float)))
        
        simulation_cell_angle = URIRef(f'{self._name}_SimulationCellAngle')
        data = self.schema.simulation_cell.angle()
        self.graph.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle))
        self.graph.add((simulation_cell_angle, RDF.type, CMSO.SimulationCellAngle))
        self.graph.add((simulation_cell_angle, CMSO.hasAngle_alpha, Literal(data[0], datatype=XSD.float)))
        self.graph.add((simulation_cell_angle, CMSO.hasAngle_beta, Literal(data[1], datatype=XSD.float)))
        self.graph.add((simulation_cell_angle, CMSO.hasAngle_gamma, Literal(data[2], datatype=XSD.float)))
        
    
    def _add_crystal_structure(self):
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
        self.graph.add((self.material, CMSO.hasStructure, crystal_structure))
        self.graph.add((crystal_structure, RDF.type, CMSO.CrystalStructure))    
        self.graph.add((crystal_structure, CMSO.hasAltName, 
            Literal(self.schema.material.crystal_structure.name(), 
                datatype=XSD.string)))
        self.crystal_structure = crystal_structure
        
    def _add_space_group(self):
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
        self.graph.add((self.crystal_structure, CMSO.hasSpaceGroup, space_group))
        self.graph.add((space_group, CMSO.hasSpaceGroupSymbol, 
            Literal(self.schema.material.crystal_structure.spacegroup_symbol(), 
                datatype=XSD.string)))
        self.graph.add((space_group, CMSO.hasSpaceGroupNumber, 
            Literal(self.schema.material.crystal_structure.spacegroup_number(), 
                datatype=XSD.integer)))
    
            
    def _add_unit_cell(self):
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
        self.graph.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell))
        self.graph.add((unit_cell, RDF.type, CMSO.UnitCell))
        self.unit_cell = unit_cell
        
        #add bravais lattice
        bv = self.schema.material.crystal_structure.unit_cell.bravais_lattice()
        if bv is not None:
            bv = URIRef(bv)
            self.graph.add((self.unit_cell, CMSO.hasBravaisLattice, bv))
        
    def _add_lattice_properties(self):
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
        data = self.schema.material.crystal_structure.unit_cell.lattice_parameter()
        lattice_parameter = URIRef(f'{self._name}_LatticeParameter')
        self.graph.add((self.unit_cell, CMSO.hasLatticeParamter, lattice_parameter))
        self.graph.add((lattice_parameter, RDF.type, CMSO.LatticeParameter))
        self.graph.add((lattice_parameter, CMSO.hasLength_x, Literal(data[0], datatype=XSD.float)))
        self.graph.add((lattice_parameter, CMSO.hasLength_y, Literal(data[1], datatype=XSD.float)))
        self.graph.add((lattice_parameter, CMSO.hasLength_z, Literal(data[2], datatype=XSD.float)))
        
        lattice_angle = URIRef(f'{self._name}_LatticeAngle')
        data = self.schema.material.crystal_structure.unit_cell.angle()
        self.graph.add((self.unit_cell, CMSO.hasAngle, lattice_angle))
        self.graph.add((lattice_angle, RDF.type, CMSO.LatticeAngle))
        self.graph.add((lattice_angle, CMSO.hasAngle_alpha, Literal(data[0], datatype=XSD.float)))
        self.graph.add((lattice_angle, CMSO.hasAngle_beta, Literal(data[1], datatype=XSD.float)))
        self.graph.add((lattice_angle, CMSO.hasAngle_gamma, Literal(data[2], datatype=XSD.float)))        


    def _save_atom_attributes(self, position_identifier, species_identifier):
        #if self.store == 'pyiron':
        #    pass
        #else:
        #    #this is the file based store system
        datadict = {
            position_identifier:{
                "value": self.schema.atom_attribute.position(),
                "label": "position", 
            },
            species_identifier:{
                "value": self.schema.atom_attribute.species(),
                "label": "species", 
            },
        }
        outfile = os.path.join(self.graph.structure_store, str(self._name).split(':')[-1])
        json_io.write_file(outfile,  datadict)
        return os.path.relpath(outfile+'.json')

    def _add_atoms(self):
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

        if "positions" in self.atoms.keys():
            position = URIRef(f'{self._name}_Position')
            self.graph.add((self.sample, CMSO.hasAttribute, position))
            self.graph.add((position, RDF.type, CMSO.AtomAttribute))
            self.graph.add((position, CMSO.hasName, Literal('Position', datatype=XSD.string)))
            self.graph.add((position, CMSO.hasIdentifier, Literal(position_identifier, datatype=XSD.string)))            
            self.graph.add((position, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))

        if "species" in self.atoms.keys():
            species = URIRef(f'{self._name}_Species')
            self.graph.add((self.sample, CMSO.hasAttribute, species))
            self.graph.add((species, RDF.type, CMSO.AtomAttribute))
            self.graph.add((species, CMSO.hasName, Literal('Species', datatype=XSD.string)))
            self.graph.add((species, CMSO.hasIdentifier, Literal(species_identifier, datatype=XSD.string)))            
            self.graph.add((species, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))

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
        if self.graph is None:
            return

        vacancy = URIRef(f'{self._name}_Vacancy')
        self.graph.add((self.material, CMSO.hasDefect, vacancy))
        self.graph.add((vacancy, RDF.type, PODO.Vacancy))
        self.graph.add((self.simulation_cell, PODO.hasVacancyConcentration, Literal(concentration, datatype=XSD.float)))
        if number is not None:
            self.graph.add((self.simulation_cell, PODO.hasNumberOfVacancies, Literal(number, datatype=XSD.integer)))

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
        if self.graph is None:
            return
                    
        if gb_dict["GBType"] is None:
            plane_defect = URIRef(f'{self._name}_GrainBoundary')
            self.graph.add((plane_defect, RDF.type, PLDO.GrainBoundary))
        
        elif gb_dict["GBType"] == "Twist":
            plane_defect = URIRef(f'{self._name}_TwistGrainBoundary')
            self.graph.add((plane_defect, RDF.type, PLDO.TwistGrainBoundary))
        
        elif gb_dict["GBType"] == "Tilt":
            plane_defect = URIRef(f'{self._name}_TiltGrainBoundary')
            self.graph.add((plane_defect, RDF.type, PLDO.TiltGrainBoundary))
        
        elif gb_dict["GBType"] == "Symmetric Tilt":
            plane_defect = URIRef(f'{self._name}_SymmetricalTiltGrainBoundary')
            self.graph.add((plane_defect, RDF.type, PLDO.SymmetricalTiltGrainBoundary))
        
        elif gb_dict["GBType"] == "Mixed":
            plane_defect = URIRef(f'{self._name}_MixedGrainBoundary')
            self.graph.add((plane_defect, RDF.type, PLDO.MixedGrainBoundary))
        
        self.graph.add((self.material, CMSO.hasDefect, plane_defect))
        self.graph.add((plane_defect, PLDO.hasSigmaValue, Literal(gb_dict["sigma"], datatype=XSD.integer)))
        self.graph.add((plane_defect, PLDO.hasGBPlane, Literal(gb_dict["GBPlane"], 
                                                             datatype=XSD.string)))
        self.graph.add((plane_defect, PLDO.hasRotationAxis, Literal(gb_dict["RotationAxis"], 
                                                             datatype=XSD.string)))
        self.graph.add((plane_defect, PLDO.hasMisorientationAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))

