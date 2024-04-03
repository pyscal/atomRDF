"""
StructureGraph is the central object in pyscal_rdf which combines all the functionality
of :py:class:`pyscal_rdf.graph.RDFGraph` along with easy structural creation routines.
"""
import numpy as np
import copy
from functools import partial, update_wrapper

from pyscal_rdf.rdfsystem import System

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

    def to_graph(self):
        if self.graph is not None:
            self._generate_name()
            self._add_sample()
            self._add_material()

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
        self.add((self.sample, CMSO.hasSpecies, chemical_species))
        self.add((chemical_species, RDF.type, CMSO.ChemicalSpecies))

        for e, r in composition.items():
            if e in element_indetifiers.keys():
                element = URIRef(element_indetifiers[e])
                self.graph.add((chemical_species, CMSO.hasElement, element))
                self.graph.add((element, RDF.type, CMSO.Element))
                self.graph.add((element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string)))
                self.graph.add((element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float)))

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
        
        self.add((self.material, CMSO.hasDefect, plane_defect))
        self.add((plane_defect, PLDO.hasSigmaValue, Literal(gb_dict["sigma"], datatype=XSD.integer)))
        self.add((plane_defect, PLDO.hasGBPlane, Literal(gb_dict["GBPlane"], 
                                                             datatype=XSD.string)))
        self.add((plane_defect, PLDO.hasRotationAxis, Literal(gb_dict["RotationAxis"], 
                                                             datatype=XSD.string)))
        self.add((plane_defect, PLDO.hasMisorientationAngle, Literal(gb_dict["MisorientationAngle"], datatype=XSD.float)))


