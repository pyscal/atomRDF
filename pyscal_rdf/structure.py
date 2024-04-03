"""
StructureGraph is the central object in pyscal_rdf which combines all the functionality
of :py:class:`pyscal_rdf.graph.RDFGraph` along with easy structural creation routines.
"""
import numpy as np
import copy
from functools import partial, update_wrapper

from pyscal_rdf.rdfsystem import System
from pyscal_rdf.graph import RDFGraph

import pyscal3.structure_creator as pcs
from pyscal3.grain_boundary import GrainBoundary
from pyscal3.atoms import AttrSetter
import pyscal_rdf.properties as prp
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict

from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("https://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")


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
    
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    if graph is not None:
        graph.add_structure_to_graph(s, names=names)
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
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = 'custom'
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    if graph is not None:
        graph.add_structure_to_graph(s, names=names)

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
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    if graph is not None:
        graph.add_structure_to_graph(s, names=names)
        gb_dict = {"GBPlane": " ".join(np.array(gb_plane).astype(str)),
                  "RotationAxis": axis,
                  "MisorientationAngle": gb.theta,
                  "GBType": gb.find_gb_character(),
                  "sigma": gb.sigma,
                  }
        graph.add_gb(gb_dict)
    return s

def _read_structure(self, filename, 
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

    s = System(filename, format=format, species=species)
    s.lattice_properties = datadict
    
    if graph is not None:
        graph.add_structure_to_graph(s, names=names)
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
        self.graph = None
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


    def __delitem__(self, val):
        if isinstance(val, int):
            val = [val]
        
        #now the graph has to be updated accordingly
        if self.graph is not None:
            #first annotate graph
            c = (len(val)/self.natoms)
            self.graph.add_vacancy(c, number=len(val))
            #now we need to re-add atoms, so at to remove
            self.graph.graph.remove((self.sample, CMSO.hasNumberOfAtoms, None))
            self.graph.graph.add((self.sample, CMSO.hasNumberOfAtoms, Literal(self.natoms-len(val), datatype=XSD.integer)))
            #revamp composition
            #remove existing chem composution
            chemical_species = self.graph.value(self.sample, CMSO.hasSpecies)
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
                    self.add((chemical_species, CMSO.hasElement, element))
                    self.add((element, RDF.type, CMSO.Element))
                    self.add((element, CMSO.hasSymbol, Literal(e, datatype=XSD.string)))
                    self.add((element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float)))
        self.delete(indices=list(val))








class StructureGraph(RDFGraph):
    def __init__(self, graph_file=None, 
        store="Memory", 
        store_file=None,
        identifier="http://default_graph",
        ontology=None,
        structure_store=None):
        
        super().__init__(graph_file=graph_file, 
            store=store, 
            store_file=store_file, 
            identifier=identifier,
            ontology=ontology,
            structure_store=structure_store)
        self._element_dict = element_dict
        self._structure_dict = structure_dict

        #add create methods
        self.create = AttrSetter()
        mapdict = {}

        #element creation routines
        mapdict["element"] = {}
        for key in self._element_dict.keys():
            mapdict["element"][key] = update_wrapper(partial(self._annotated_make_crystal,
                self._element_dict[key]['structure'],
                lattice_constant=self._element_dict[key]['lattice_constant'],
                element = key), pcs.make_crystal)
        
        #lattice creation routines
        mapdict["lattice"] = {}
        for key in self._structure_dict.keys():
            mapdict["lattice"][key] = update_wrapper(partial(self._annotated_make_crystal, 
                key), 
                _make_crystal)

        #custom creation routines
        mapdict["lattice"]["custom"] = self._annotated_make_general_lattice

        #create defects
        mapdict["defect"] = {}
        mapdict["defect"]["grain_boundary"] = self._annotated_make_grain_boundary

        self.create._add_attribute(mapdict)



    def _annotated_make_crystal(self, structure, 
                lattice_constant = 1.00, 
                repetitions = None, 
                ca_ratio = 1.633, 
                noise = 0, 
                element=None,
                primitive=False,
                add_to_graph=True, 
                names=False):

        sys = System(source=_make_crystal(structure, 
                lattice_constant = lattice_constant, 
                repetitions = repetitions, 
                ca_ratio = ca_ratio, 
                noise = noise, 
                element=element,
                primitive=primitive))
        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)
        return sys

    def _annotated_make_general_lattice(self, positions,
            types, 
            box,
            lattice_constant = 1.00, 
            repetitions = None, 
            noise = 0,
            element=None,
            add_to_graph=True,
            names=False):
        
        sys = System(source=_make_general_lattice(positions,
            types, 
            box,
            lattice_constant = lattice_constant, 
            repetitions = repetitions, 
            noise = noise,
            element=element))

        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)
        return sys

    def _annotated_make_grain_boundary(self, axis, 
            sigma, gb_plane,
            structure = None,
            element = None, 
            lattice_constant = 1,
            repetitions = (1,1,1),
            overlap = 0.0,
            add_to_graph = True,
            names = False):
        """
        Create a grain boundary structure and return it as a System object.

        Parameters
        ----------
        axis: list of ints of length 3
            The grain boundary axis

        sigma: int
            sigma value of the grain boundary
        
        gb_plane: list of ints of length 3
            The grain boundary plane

        structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
            type of the crystal structure

        element : string, optional
            The chemical element

        lattice_constant : float, optional
            lattice constant of the crystal structure, default 1

        repetitions : list of ints of len 3, optional
            of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
            default `[1, 1, 1]`.

        overlap: float, optional
            overlap between the two grains

        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs            

        Returns
        -------
        System: pyscal System

        """

        sys, gb = _make_grain_boundary(axis, 
            sigma, gb_plane,
            structure = structure,
            element = element, 
            lattice_constant = lattice_constant,
            repetitions = repetitions,
            overlap = overlap)
        sys = System(source=sys)

        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)

        gb_dict = {"GBPlane": " ".join(np.array(gb_plane).astype(str)),
                  "RotationAxis": axis,
                  "MisorientationAngle": gb.theta,
                  "GBType": gb.find_gb_character(),
                  "sigma": gb.sigma,
                  }
        self.add_gb(gb_dict)

        return sys


    
    def read_structure(self, filename, format="lammps-dump",
                      add_to_graph=True, names=False,
                      species=None,
                      lattice=None,
                      lattice_constant=None,
                      basis_box=None,
                      basis_positions=None,
                      ):
        """
        Read an input file and return it as a System object.

        Parameters
        ----------
        filename: string
            name of the input file

        format: string
            format of the input file

        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs

        Returns
        -------
        System: pyscal System
        system will be populated with given atoms and simulation box

        """
        #prepare the dict for storing extra info; if lattice is provided, extract it
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

        sys = System(filename, format=format, species=species)
        sys.lattice_properties = datadict
        
        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)
            #sys.sample = self.sample
            #sys._atom_ids = copy.copy(self._atom_ids)
            #sys.graph = self
        return sys            

