import numpy as np
from functools import partial, update_wrapper

import pyscal3.core as pc
from pyscal3.atoms import AttrSetter
import pyscal_rdf.properties as prp

from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("https://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")

class System(pc.System):
    def __init__(self, filename = None, 
            format = "lammps-dump", 
            compressed = False, 
            customkeys = None,
            species = None,
            source=None):
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