"""
This module provides functions for creating and manipulating atomic structures. It includes functions for creating crystals,
general lattices, dislocations, grain boundaries, and reading structures from files. The module also provides functionality for
adding interstitial impurities, substituting atoms, deleting atoms, and adding vacancies.
The structures can be converted to RDF graphs using the atomrdf library.
The main object in this module is the System class, which extends the functionality of the pyscal3.core.System class and provides additional methods for working with atomic structures.
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
import warnings
from ase.io import write
from ase.build import cut

import pyscal3.structure_creator as pcs
from pyscal3.grain_boundary import GrainBoundary
from pyscal3.atoms import AttrSetter, Atoms
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict
from pyscal3.formats.ase import convert_snap

import pyscal3.operations.input as inputmethods
import pyscal3.operations.serialize as serialize
import pyscal3.operations.visualize as visualize
import pyscal3.operations.operations as operations

import atomrdf.json_io as json_io
import atomrdf.properties as prp
from atomrdf.sample import Property
import atomrdf.io as aio

from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.namespace import CMSO, LDO, PLDO, PODO, CDCO, UNSAFEASMO, UNSAFECMSO, PROV, Literal, ASMO

# read element data file
file_location = os.path.dirname(__file__).split("/")
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__), "data/element.yml")
with open(file_location, "r") as fin:
    element_indetifiers = yaml.safe_load(fin)
    
class System(pc.System):
    def __init__(
        self,
        filename=None,
        format="lammps-dump",
        compressed=False,
        customkeys=None,
        species=None,
        source=None,
        graph=None,
        names=False,
        warn_read_in=True,
    ):

        if (filename is not None) and warn_read_in:
            warnings.warn(
                "To provide additional information, use the read method"
            )

        super().__init__(
            filename=filename,
            format=format,
            compressed=compressed,
            customkeys=customkeys,
            species=species,
        )

        # this is the sample which will be stored
        self.sample = None
        self.label = None
        # the graph object should also be attached
        # for post-processing of structures
        self.graph = graph
        self.names = names
        self._material = None
        self._name = None
        self._atom_ids = None
        if source is not None:
            self.__dict__.update(source.__dict__)

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
                "repetitions": partial(prp.get_repetitions, self),
            },
            "atom_attribute": {
                "position": partial(prp.get_position, self),
                "species": partial(prp.get_species, self),
            },
        }

        self.schema._add_attribute(mapdict)


    @property
    def material(self):
        if self._material is None:
            self._material = self.graph.value(self.sample, CMSO.hasMaterial)
        return self._material

    @material.setter
    def material(self, value):
        self._material = value

    def duplicate(self, only_essential=False):
        new_system = System()
        if only_essential:
            n_dict = {'positions': copy.deepcopy(self.atoms.positions),
                    'species': copy.deepcopy(self.atoms.species),
                    'types': copy.deepcopy(self.atoms.types),}
        else:
            n_dict = {key: copy.deepcopy(val)[:self.natoms] for key, val in self.atoms.items()}
            new_system.label = self.label
            new_system._name = self._name
        
        atoms = Atoms()
        atoms.from_dict(n_dict)
        atoms._lattice = self.atoms._lattice
        atoms._lattice_constant = self.atoms._lattice_constant
        new_system._structure_dict = copy.deepcopy(self._structure_dict)
        new_system.box = self.box
        new_system.atoms = atoms
        new_system.graph = self.graph
        new_system.sample = None
        return new_system

        
    def add_property_mappings(self, output_property, mapping_quantity=None):
        if self.graph is None:
            return
        if not isinstance(output_property, Property):
            return
        
        #if the property is a directly calculated value
        parent_samples = list([x[0] for x in self.graph.triples((None, ASMO.hasCalculatedProperty, output_property._parent))])
        if len(parent_samples)>0:
            for parent_sample in parent_samples:
                self.graph.add((self.sample, PROV.wasDerivedFrom, parent_sample))
        else:
            #this is quantity that is derived -> for example volume/3 -> it has only sample parent, but no direct connection
            if output_property._sample_parent is not None:
                self.graph.add((self.sample, PROV.wasDerivedFrom, output_property._sample_parent))

        if mapping_quantity=='lattice_constant':
            #add lattice constant mapping
            material = self.graph.value(self.sample, CMSO.hasMaterial)
            crystal_structure = self.graph.value(material, CMSO.hasStructure)
            unit_cell = self.graph.value(crystal_structure, CMSO.hasUnitCell)
            lattice_parameter = self.graph.value(unit_cell, CMSO.hasLatticeParameter)

            #also get activity
            activity = self.graph.value(output_property._parent, ASMO.wasCalculatedBy)
            self.graph.add((lattice_parameter, ASMO.wasCalculatedBy, activity))

    def add_vacancy(self, concentration, number=None):
        """
        Add Vacancy details which will be annotated by PODO

        Parameters
        ----------
        concentration: float
            vacancy concentration, value should be between 0-1

        number: int
            Number of atoms that were deleted, optional

        Returns
        -------
        None
        """
        if self.graph is None:
            return

        vacancy = self.graph.create_node(f"{self._name}_Vacancy", PODO.Vacancy)
        self.graph.add((self.material, CDCO.hasCrystallographicDefect, vacancy))
        self.graph.add(
            (
                self.sample,
                PODO.hasVacancyConcentration,
                Literal(concentration, datatype=XSD.float),
            )
        )
        if number is not None:
            self.graph.add(
                (
                    self.sample,
                    PODO.hasNumberOfVacancies,
                    Literal(number, datatype=XSD.integer),
                )
            )

    def update_system_for_vacancy_creation(self, vacancy_no, actual_natoms, original_sample):
        if self.graph is None:
            return

        # now we need to re-add atoms, so at to remove
        self.graph.remove((self.sample, CMSO.hasNumberOfAtoms, None))
        self.graph.add(
            (
                self.sample,
                CMSO.hasNumberOfAtoms,
                Literal(actual_natoms - vacancy_no, datatype=XSD.integer),
            )
        )

        chemical_species = self.graph.value(self.sample, CMSO.hasSpecies)
        # start by cleanly removing elements
        for s in self.graph.triples((chemical_species, CMSO.hasElement, None)):
            element = s[2]
            self.graph.remove((element, None, None))
        self.graph.remove((chemical_species, None, None))
        self.graph.remove((self.sample, CMSO.hasSpecies, None))

        # now recalculate and add it again
        composition = self.schema.material.element_ratio()
        valid = False
        for e, r in composition.items():
            if e in element_indetifiers.keys():
                valid = True
                break

        if valid:
            chemical_species = self.graph.create_node(
                f"{self._name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = self.graph.create_node(
                        element_indetifiers[e], CMSO.ChemicalElement
                    )
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add(
                        (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                    )
                    self.graph.add(
                        (
                            element,
                            CMSO.hasElementRatio,
                            Literal(r, datatype=XSD.float),
                        )
                    )

        # we also have to read in file and clean it up
        filepath = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasPath
        ).toPython()
        position_identifier = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasIdentifier
        ).toPython()
        species_identifier = self.graph.value(
            URIRef(f"{self.sample}_Species"), CMSO.hasIdentifier
        ).toPython()

        # clean up items
        datadict = {
            position_identifier: {
                "value": self.schema.atom_attribute.position(),
                "label": "position",
            },
            species_identifier: {
                "value": self.schema.atom_attribute.species(),
                "label": "species",
            },
        }
        outfile = os.path.join(
            self.graph.structure_store, str(self._name).split(":")[-1]
        )
        json_io.write_file(outfile, datadict)

        #write mapping for the operation
        if original_sample != self.sample:
            activity = self.graph.create_node(f"structuremanipulation:{uuid.uuid4()}", ASMO.DeleteAtom)
            self.graph.add((self.sample, PROV.wasDerivedFrom, original_sample))
            self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

    def add_substitutional_impurities(self, 
                conc_of_impurities, 
                no_of_impurities=None):
        defect = self.graph.create_node(f"{self._name}_SubstitutionalImpurity", PODO.SubstitutionalImpurity)
        self.graph.add((self.material, CDCO.hasCrystallographicDefect, defect))
        self.graph.add((self.sample, PODO.hasImpurityConcentration, Literal(conc_of_impurities, datatype=XSD.float)))
        if no_of_impurities is not None:
            self.graph.add((self.sample, PODO.hasNumberOfImpurityAtoms, Literal(no_of_impurities, datatype=XSD.integer)))

    def update_system_for_substitutional_impurity(self, 
                impurity_no, 
                actual_natoms,
                original_sample):
        # operate on the graph
        if self.graph is None:
            return
        
        chemical_species = self.graph.value(self.sample, CMSO.hasSpecies)
        # start by cleanly removing elements
        for s in self.graph.triples((chemical_species, CMSO.hasElement, None)):
            element = s[2]
            self.graph.remove((element, None, None))
        self.graph.remove((chemical_species, None, None))
        self.graph.remove((self.sample, CMSO.hasSpecies, None))

        composition = self.schema.material.element_ratio()
        valid = False
        for e, r in composition.items():
            if e in element_indetifiers.keys():
                valid = True
                break

        if valid:
            chemical_species = self.graph.create_node(
                f"{self._name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = self.graph.create_node(
                        element_indetifiers[e], CMSO.ChemicalElement
                    )
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add(
                        (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                    )
                    self.graph.add(
                        (
                            element,
                            CMSO.hasElementRatio,
                            Literal(r, datatype=XSD.float),
                        )
                    )

        # we also have to read in file and clean it up
        filepath = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasPath
        ).toPython()
        position_identifier = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasIdentifier
        ).toPython()
        species_identifier = self.graph.value(
            URIRef(f"{self.sample}_Species"), CMSO.hasIdentifier
        ).toPython()

        # clean up items
        datadict = {
            position_identifier: {
                "value": self.schema.atom_attribute.position(),
                "label": "position",
            },
            species_identifier: {
                "value": self.schema.atom_attribute.species(),
                "label": "species",
            },
        }
        outfile = os.path.join(
            self.graph.structure_store, str(self._name).split(":")[-1]
        )
        json_io.write_file(outfile, datadict)

        #write mapping for the operation
        if original_sample != self.sample:
            activity = self.graph.create_node(f"activity:{uuid.uuid4()}", ASMO.SubstituteAtom)
            self.graph.add((self.sample, PROV.wasDerivedFrom, original_sample))
            self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

    def add_interstitial_impurities(self, 
            conc_of_impurities, 
            no_of_impurities=None, 
            label=None):
        if label is not None:
            defect = self.graph.create_node(f"{self._name}_InterstitialImpurity", PODO.InterstitialImpurity, label=label)
        else:
            defect = self.graph.create_node(f"{self._name}_InterstitialImpurity", PODO.InterstitialImpurity)
        self.graph.add((self.material, CDCO.hasCrystallographicDefect, defect))
        self.graph.add((self.sample, PODO.hasImpurityConcentration, Literal(conc_of_impurities, datatype=XSD.float)))
        if no_of_impurities is not None:
            self.graph.add((self.sample, PODO.hasNumberOfImpurityAtoms, Literal(no_of_impurities, datatype=XSD.integer)))

    def update_system_for_interstitial_impurity(self,
        original_sample):
        if self.graph is None:
            return
        self.graph.remove((self.sample, CMSO.hasNumberOfAtoms, None))
        self.graph.add(
            (
                self.sample,
                CMSO.hasNumberOfAtoms,
                Literal(self.natoms, datatype=XSD.integer),
            )
        )
        # revamp composition
        # remove existing chem composution
        chemical_species = self.graph.value(self.sample, CMSO.hasSpecies)
        # start by cleanly removing elements
        for s in self.graph.triples((chemical_species, CMSO.hasElement, None)):
            element = s[2]
            self.graph.remove((element, None, None))
        self.graph.remove((chemical_species, None, None))
        self.graph.remove((self.sample, CMSO.hasSpecies, None))

        composition = self.schema.material.element_ratio()
        valid = False
        for e, r in composition.items():
            if e in element_indetifiers.keys():
                valid = True
                break

        if valid:
            chemical_species = self.graph.create_node(
                f"{self._name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = self.graph.create_node(
                        element_indetifiers[e], CMSO.ChemicalElement
                    )
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add(
                        (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                    )
                    self.graph.add(
                        (
                            element,
                            CMSO.hasElementRatio,
                            Literal(r, datatype=XSD.float),
                        )
                    )

        # we also have to read in file and clean it up
        filepath = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasPath
        ).toPython()
        position_identifier = self.graph.value(
            URIRef(f"{self.sample}_Position"), CMSO.hasIdentifier
        ).toPython()
        species_identifier = self.graph.value(
            URIRef(f"{self.sample}_Species"), CMSO.hasIdentifier
        ).toPython()

        # clean up items
        datadict = {
            position_identifier: {
                "value": self.schema.atom_attribute.position(),
                "label": "position",
            },
            species_identifier: {
                "value": self.schema.atom_attribute.species(),
                "label": "species",
            },
        }
        outfile = os.path.join(
            self.graph.structure_store, str(self._name).split(":")[-1]
        )
        json_io.write_file(outfile, datadict)

        #write mapping for the operation
        if original_sample != self.sample:
            activity = self.graph.create_node(f"activity:{uuid.uuid4()}", ASMO.AddAtom)
            self.graph.add((self.sample, PROV.wasDerivedFrom, self.sample))
            self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

    def to_graph(self):
        """
        Converts the structure object to a graph representation.

        Returns
        -------
        None
        """
        if self.graph is None:
            return

        self._generate_name()
        self._add_sample()
        self._add_material()
        self._add_chemical_composition()
        self._add_simulation_cell()
        self._add_simulation_cell_properties()
        self._add_crystal_structure()
        self._add_atoms()

    def _generate_name(self, name_index=None):
        if self.names:
            if name_index is None:
                name_index = self.graph.n_samples + 1
            self._name = f"sample:{name_index}"
        else:
            self._name = f"sample:{str(uuid.uuid4())}"

    def _add_sample(self):
        sample = self.graph.create_node(self._name, CMSO.AtomicScaleSample, label=self.label)
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
        material = self.graph.create_node(
            f"{self._name}_Material", CMSO.CrystallineMaterial
        )
        self.graph.add((self.sample, CMSO.hasMaterial, material))
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
        valid = False
        for e, r in composition.items():
            if e in element_indetifiers.keys():
                valid = True
                break

        if valid:
            chemical_species = self.graph.create_node(
                f"{self._name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = self.graph.create_node(
                        element_indetifiers[e], CMSO.ChemicalElement
                    )
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add(
                        (
                            element,
                            CMSO.hasChemicalSymbol,
                            Literal(e, datatype=XSD.string),
                        )
                    )
                    self.graph.add(
                        (element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float))
                    )

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

        simulation_cell = self.graph.create_node(
            f"{self._name}_SimulationCell", CMSO.SimulationCell
        )
        self.graph.add((self.sample, CMSO.hasSimulationCell, simulation_cell))
        volume = self.graph.create_node(
            f"{self._name}_Volume", ASMO.Volume, label="SimulationCellVolume"
        )
        self.graph.add((simulation_cell, CMSO.hasVolume, volume))
        self.graph.add(
            (
                volume,
                ASMO.hasValue,
                Literal(
                    np.round(self.schema.simulation_cell.volume(), decimals=2),
                    datatype=XSD.float,
                ),
            )
        )
        self.graph.add(
            (
                volume,
                ASMO.hasUnit,
                URIRef(f"http://qudt.org/vocab/unit/ANGSTROM3"),
            )
        )
        self.graph.add(
            (
                self.sample,
                CMSO.hasNumberOfAtoms,
                Literal(
                    self.schema.simulation_cell.number_of_atoms(), datatype=XSD.integer
                ),
            )
        )
        
        repetitions = self.schema.simulation_cell.repetitions()
        self.graph.add((simulation_cell, CMSO.hasRepetition_x, Literal(repetitions[0], datatype=XSD.integer)))
        self.graph.add((simulation_cell, CMSO.hasRepetition_y, Literal(repetitions[1], datatype=XSD.integer)))
        self.graph.add((simulation_cell, CMSO.hasRepetition_z, Literal(repetitions[2], datatype=XSD.integer)))
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
        simulation_cell_length = self.graph.create_node(
            f"{self._name}_SimulationCellLength", CMSO.SimulationCellLength
        )
        self.graph.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length))
        data = self.schema.simulation_cell.length()
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_x,
                Literal(data[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_y,
                Literal(data[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_z,
                Literal(data[2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_01 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_1", CMSO.SimulationCellVector
        )
        data = self.schema.simulation_cell.vector()
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01)
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_x,
                Literal(data[0][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_y,
                Literal(data[0][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_z,
                Literal(data[0][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_02 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_2", CMSO.SimulationCellVector
        )
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02)
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_x,
                Literal(data[1][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_y,
                Literal(data[1][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_z,
                Literal(data[1][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_03 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_3", CMSO.SimulationCellVector
        )
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03)
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_x,
                Literal(data[2][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_y,
                Literal(data[2][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_z,
                Literal(data[2][2], datatype=XSD.float),
            )
        )

        simulation_cell_angle = self.graph.create_node(
            f"{self._name}_SimulationCellAngle", CMSO.SimulationCellAngle
        )
        data = self.schema.simulation_cell.angle()
        self.graph.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle))
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_alpha,
                Literal(data[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_beta,
                Literal(data[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_gamma,
                Literal(data[2], datatype=XSD.float),
            )
        )

    def _add_crystal_structure(self, targets=None):
        """
        Add a CMSO Crystal Structure

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        if targets is None:
            targets = [
                self.schema.material.crystal_structure.name(),
                self.schema.material.crystal_structure.spacegroup_symbol(),
                self.schema.material.crystal_structure.spacegroup_number(),
                self.schema.material.crystal_structure.unit_cell.bravais_lattice(),
                self.schema.material.crystal_structure.unit_cell.lattice_parameter(),
                self.schema.material.crystal_structure.unit_cell.angle(),
            ]

        #fix for lattice angle of HCP
        if targets[0] == 'hcp':
            targets[5] = [90.0, 90.0, 120.0]

        valid = self.graph._is_valid(targets)

        if valid:
            crystal_structure = self.graph.create_node(
                f"{self._name}_CrystalStructure", CMSO.CrystalStructure
            )
            self.graph.add((self.material, CMSO.hasStructure, crystal_structure))
            self.graph.add(
                (
                    crystal_structure,
                    CMSO.hasAltName,
                    Literal(targets[0], datatype=XSD.string),
                )
            )
            self.crystal_structure = crystal_structure

            if targets[1] is not None:
                self._add_space_group(targets[1], targets[2])

            # now see if unit cell needs to be added
            valid = self.graph._is_valid(targets[3:])
            if valid:
                self._add_unit_cell()
                if targets[3] is not None:
                    self._add_bravais_lattice(targets[3])
                if targets[4] is not None:
                    self._add_lattice_properties(targets[4], targets[5])

    def _add_space_group(self, spacegroup_symbol, spacegroup_number):
        """
        Add a CMSO Space Group

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        #space_group = URIRef(f"{self._name}_SpaceGroup")
        #self.graph.add((self.crystal_structure, CMSO.hasSpaceGroup, space_group))
        self.graph.add(
            (
                self.crystal_structure,
                CMSO.hasSpaceGroupSymbol,
                Literal(spacegroup_symbol, datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                self.crystal_structure,
                CMSO.hasSpaceGroupNumber,
                Literal(spacegroup_number, datatype=XSD.integer),
            )
        )

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

        unit_cell = self.graph.create_node(f"{self._name}_UnitCell", CMSO.UnitCell)
        self.graph.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell))
        self.unit_cell = unit_cell

    def _add_bravais_lattice(self, bv):
        """
        Add a Bravais lattice to the unit cell.

        Parameters:
            bv (str): The URI of the Bravais lattice.

        Returns:
            None
        """
        bv = URIRef(bv)
        self.graph.add(
            (
                self.unit_cell,
                CMSO.hasBravaisLattice,
                bv,
            )
        )

    def _add_lattice_properties(self, lattice_parameter_value, lattice_angle_value):
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
        lattice_parameter = self.graph.create_node(
            f"{self._name}_LatticeParameter", CMSO.LatticeParameter
        )
        self.graph.add(
            (
                self.unit_cell,
                CMSO.hasLength_x,
                Literal(lattice_parameter_value[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                self.unit_cell,
                CMSO.hasLength_y,
                Literal(lattice_parameter_value[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                self.unit_cell,
                CMSO.hasLength_z,
                Literal(lattice_parameter_value[2], datatype=XSD.float),
            )
        )

        lattice_angle = self.graph.create_node(
            f"{self._name}_LatticeAngle", CMSO.LatticeAngle
        )
        self.graph.add((self.unit_cell, CMSO.hasAngle, lattice_angle))
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_alpha,
                Literal(lattice_angle_value[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_beta,
                Literal(lattice_angle_value[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_gamma,
                Literal(lattice_angle_value[2], datatype=XSD.float),
            )
        )

    def _save_atom_attributes(self, position_identifier, species_identifier):
        """
        Save the atom attributes to a file.

        Parameters
        ----------
        position_identifier : str
            The identifier for the position attribute.
        species_identifier : str
            The identifier for the species attribute.

        Returns
        -------
        str
            The relative path to the saved file.

        Notes
        -----
        This method saves the atom attributes to a file in the file-based store system.
        The attributes are stored in a dictionary with the position identifier and species identifier as keys.
        The dictionary is then written to a JSON file using the `json_io.write_file` function.
        The file is saved in the structure store directory with the name of the structure as the filename.
        The method returns the relative path to the saved file.
        """
        datadict = {
            position_identifier: {
                "value": self.schema.atom_attribute.position(),
                "label": "position",
            },
            species_identifier: {
                "value": self.schema.atom_attribute.species(),
                "label": "species",
            },
        }
        outfile = os.path.join(
            self.graph.structure_store, str(self._name).split(":")[-1]
        )
        json_io.write_file(outfile, datadict)
        return os.path.relpath(outfile + ".json")

    def _add_atoms(self):
        """
        Add Atoms including their species and positions

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Note that for the moment, we will dump the structures in a given folder,
        maybe this could be input from the Job class directly
        """
        # now we write out file
        position_identifier = str(uuid.uuid4())
        species_identifier = str(uuid.uuid4())

        outfile = self._save_atom_attributes(position_identifier, species_identifier)

        if "positions" in self.atoms.keys():
            position = self.graph.create_node(
                f"{self._name}_Position", CMSO.AtomAttribute
            )
            self.graph.add(
                (
                    self.sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    position,
                )
            )
            self.graph.add(
                (position, CMSO.hasName, Literal("Position", datatype=XSD.string))
            )
            self.graph.add(
                (
                    position,
                    CMSO.hasIdentifier,
                    Literal(position_identifier, datatype=XSD.string),
                )
            )
            self.graph.add(
                (position, CMSO.hasPath, Literal(outfile, datatype=XSD.string))
            )

        if "species" in self.atoms.keys():
            species = self.graph.create_node(
                f"{self._name}_Species", CMSO.AtomAttribute
            )
            self.graph.add(
                (
                    self.sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    species,
                )
            )
            self.graph.add(
                (species, CMSO.hasName, Literal("Species", datatype=XSD.string))
            )
            self.graph.add(
                (
                    species,
                    CMSO.hasIdentifier,
                    Literal(species_identifier, datatype=XSD.string),
                )
            )
            self.graph.add(
                (species, CMSO.hasPath, Literal(outfile, datatype=XSD.string))
            )

        # if "velocities" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Velocity'
        #    velocity = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, velocity))
        #    self.add((velocity, RDF.type, CMSO.AtomAttribute))
        #    self.add((velocity, CMSO.hasName, Literal('Velocity', data_type=XSD.string)))
        #    velocity_identifier = uuid.uuid4()
        #    self.add((velocity, CMSO.hasIdentifier, Literal(velocity_identifier, datatype=XSD.string)))

        # if "forces" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Force'
        #    force = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, force))
        #    self.add((force, RDF.type, CMSO.AtomAttribute))
        #    self.add((force, CMSO.hasName, Literal('Force', data_type=XSD.string)))
        #    force_identifier = uuid.uuid4()
        #    self.add((force, CMSO.hasIdentifier, Literal(force_identifier, datatype=XSD.string)))

    def add_dislocation(self, disl_dict):
        if self.graph is None:
            return
        
        #find what kind of disl is present
        angle_deg = disl_dict['DislocationCharacter']
        if (np.abs(angle_deg-0) < 1E-3) or (np.abs(angle_deg-180) < 1E-3) or (np.abs(angle_deg-360) < 1E-3):
            disl_type = LDO.ScrewDislocation
            disl_name = "ScrewDislocation"
        elif (np.abs(angle_deg-90) < 1E-3) or (np.abs(angle_deg-270) < 1E-3):
            disl_type = LDO.EdgeDislocation
            disl_name = "EdgeDislocation"
        else:
            disl_type = LDO.MixedDislocation
            disl_name = "MixedDislocation"

        line_defect = self.graph.create_node(f"{self._name}_Dislocation", disl_type)
        self.graph.add((self.material, CDCO.hasCrystallographicDefect, line_defect))

        line_direction = self.graph.create_node(f"{self._name}_DislocationLineDirection", LDO.LineDirection)
        self.graph.add((line_direction, CMSO.hasComponent_x, Literal(disl_dict['DislocationLine'][0], datatype=XSD.float)))
        self.graph.add((line_direction, CMSO.hasComponent_y, Literal(disl_dict['DislocationLine'][1], datatype=XSD.float)))
        self.graph.add((line_direction, CMSO.hasComponent_z, Literal(disl_dict['DislocationLine'][2], datatype=XSD.float)))
        self.graph.add((line_defect, LDO.hasLineDirection, line_direction))                

        burgers_vector = self.graph.create_node(f"{self._name}_DislocationBurgersVector", LDO.BurgersVector)
        self.graph.add((burgers_vector, CMSO.hasComponent_x, Literal(disl_dict['BurgersVector'][0], datatype=XSD.float)))
        self.graph.add((burgers_vector, CMSO.hasComponent_y, Literal(disl_dict['BurgersVector'][1], datatype=XSD.float)))
        self.graph.add((burgers_vector, CMSO.hasComponent_z, Literal(disl_dict['BurgersVector'][2], datatype=XSD.float)))
        self.graph.add((line_defect, LDO.hasBurgersVector, burgers_vector))

        if disl_name == "MixedDislocation":
            self.graph.add((line_defect, LDO.hasCharacterAngle, Literal(angle_deg, datatype=XSD.float)))

        slip_direction = self.graph.create_node(f"{self._name}_DislocationSlipDirection", LDO.SlipDirection)
        self.graph.add((slip_direction, CMSO.hasComponent_x, Literal(disl_dict['SlipDirection'][0], datatype=XSD.float)))
        self.graph.add((slip_direction, CMSO.hasComponent_y, Literal(disl_dict['SlipDirection'][1], datatype=XSD.float)))
        self.graph.add((slip_direction, CMSO.hasComponent_z, Literal(disl_dict['SlipDirection'][2], datatype=XSD.float)))
        
        slip_plane = self.graph.create_node(f"{self._name}_DislocationSlipPlane", LDO.SlipPlane)
        normal_vector = self.graph.create_node(f"{self._name}_DislocationNormalVector", LDO.NormalVector)
        self.graph.add((normal_vector, CMSO.hasComponent_x, Literal(disl_dict['SlipPlane'][0], datatype=XSD.float)))
        self.graph.add((normal_vector, CMSO.hasComponent_y, Literal(disl_dict['SlipPlane'][1], datatype=XSD.float)))
        self.graph.add((normal_vector, CMSO.hasComponent_z, Literal(disl_dict['SlipPlane'][2], datatype=XSD.float)))
        self.graph.add((slip_plane, LDO.hasNormalVector, normal_vector))

        slip_system = self.graph.create_node(f"{self._name}_DislocationSlipSystem", LDO.SlipSystem)
        self.graph.add((slip_direction, LDO.belongsToSystem, slip_system))
        self.graph.add((slip_plane, LDO.belongsToSystem, slip_system))
        self.graph.add((line_defect, LDO.movesOn, slip_system))


    def add_stacking_fault(self, sf_dict):
        if self.graph is None:
            return
        plane = " ".join(np.array(sf_dict["plane"]).astype(str))
        displ = " ".join(np.array(sf_dict["displacement"]).astype(str))
        sf = self.graph.create_node(f"{self._name}_StackingFault", PLDO.StackingFault)
        self.graph.add((self.material, CDCO.hasCrystallographicDefect, sf))
        self.graph.add((sf, PLDO.hasSFplane, Literal(plane, datatype=XSD.string)))
        self.graph.add((sf, PLDO.hasDisplacementVector, Literal(displ, datatype=XSD.string)))

    def add_gb(self, gb_dict):
        """
        Add GB details which will be annotated using PLDO

        Parameters
        ----------
        gb_dict : dict
            A dictionary containing details about the grain boundary.
            It should have the following keys:
            - "GBType" (str): The type of grain boundary. Possible values are "Twist", "Tilt", "Symmetric Tilt", and "Mixed".
            - "sigma" (int): The sigma value of the grain boundary.
            - "GBPlane" (str): The plane of the grain boundary.
            - "RotationAxis" (list): The rotation axis of the grain boundary.
            - "MisorientationAngle" (float): The misorientation angle of the grain boundary.

        Returns
        -------
        None

        Notes
        -----
        This method adds grain boundary details to the structure and annotates it using PLDO ontology.
        The grain boundary type, sigma value, GB plane, rotation axis, and misorientation angle are stored as attributes of the grain boundary node in the graph.
        """
        # mark that the structure has a defect
        if self.graph is None:
            return

        if gb_dict["GBType"] is None:
            plane_defect = self.graph.create_node(f"{self._name}_GrainBoundary", PLDO.GrainBoundary)

        elif gb_dict["GBType"] == "Twist":
            plane_defect = self.graph.create_node(
                f"{self._name}_TwistGrainBoundary", PLDO.TwistGrainBoundary
            )

        elif gb_dict["GBType"] == "Tilt":
            plane_defect = self.graph.create_node(
                f"{self._name}_TiltGrainBoundary", PLDO.TiltGrainBoundary
            )

        elif gb_dict["GBType"] == "Symmetric Tilt":
            plane_defect = self.graph.create_node(
                f"{self._name}_SymmetricalTiltGrainBoundary",
                PLDO.SymmetricalTiltGrainBoundary,
            )

        elif gb_dict["GBType"] == "Mixed":
            plane_defect = self.graph.create_node(
                f"{self._name}_MixedGrainBoundary", PLDO.MixedGrainBoundary
            )

        self.graph.add((self.material, CDCO.hasCrystallographicDefect, plane_defect))
        self.graph.add(
            (
                plane_defect,
                PLDO.hasSigmaValue,
                Literal(gb_dict["sigma"], datatype=XSD.integer),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasGBplane,
                Literal(gb_dict["GBPlane"], datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasRotationAxis,
                Literal(gb_dict["RotationAxis"], datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasMisorientationAngle,
                Literal(gb_dict["MisorientationAngle"], datatype=XSD.float),
            )
        )

    def add_rotation_triples(self, rotation_vectors, child_sample_id):
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Rotation)
        self.graph.add((child_sample_id, PROV.wasGeneratedBy, activity))
        self.graph.add((child_sample_id, PROV.wasDerivedFrom, self.sample))

        rot_vector_01 = self.graph.create_node(f"{activity_id}_RotationVector_1", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_01))
        self.graph.add((rot_vector_01, CMSO.hasComponent_x, Literal(rotation_vectors[0][0], datatype=XSD.float),))
        self.graph.add((rot_vector_01, CMSO.hasComponent_y, Literal(rotation_vectors[0][1], datatype=XSD.float),))
        self.graph.add((rot_vector_01, CMSO.hasComponent_z, Literal(rotation_vectors[0][2], datatype=XSD.float),))

        rot_vector_02 = self.graph.create_node(f"{activity_id}_RotationVector_2", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_02))
        self.graph.add((rot_vector_02, CMSO.hasComponent_x, Literal(rotation_vectors[1][0], datatype=XSD.float),))
        self.graph.add((rot_vector_02, CMSO.hasComponent_y, Literal(rotation_vectors[1][1], datatype=XSD.float),))
        self.graph.add((rot_vector_02, CMSO.hasComponent_z, Literal(rotation_vectors[1][2], datatype=XSD.float),))

        rot_vector_03 = self.graph.create_node(f"{activity_id}_RotationVector_3", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_03))
        self.graph.add((rot_vector_03, CMSO.hasComponent_x, Literal(rotation_vectors[2][0], datatype=XSD.float),))
        self.graph.add((rot_vector_03, CMSO.hasComponent_y, Literal(rotation_vectors[2][1], datatype=XSD.float),))
        self.graph.add((rot_vector_03, CMSO.hasComponent_z, Literal(rotation_vectors[2][2], datatype=XSD.float),))

    def _select_by_plane(self, plane, distance, reverse_orientation=False):
        plane_norm = np.linalg.norm(plane)
        selection = []
        for pos in self.atoms.positions:
            dist = np.dot(plane, pos)/plane_norm
            
            if dist < distance:
                selection.append(True)
            else:
                selection.append(False)
        if reverse_orientation:
            selection = np.invert(selection)
        return selection             

    def select_by_plane(self, plane, distance, reverse_orientation=False):
        selection = self._select_by_plane(plane, distance, 
                        reverse_orientation=reverse_orientation)
        self.apply_selection(condition=selection)
        
    
    def add_translation_triples(self, translation_vector, plane, distance, original_sample):
        if self.graph is None:
            return
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Translation)
        self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

        #now add specifics
        #shear is a vector
        t_vector = self.graph.create_node(f"{activity_id}_TranslationVector", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, t_vector))
        self.graph.add((t_vector, CMSO.hasComponent_x, Literal(translation_vector[0], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_y, Literal(translation_vector[1], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_z, Literal(translation_vector[2], datatype=XSD.float),))

        if self.sample != original_sample:
            sys.graph.add((sys.sample, PROV.wasDerivedFrom, original_sample))
    
    def add_shear_triples(self, translation_vector, plane, distance, original_sample):
        if self.graph is None:
            return
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Shear)
        self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

        #now add specifics
        #shear is a vector
        t_vector = self.graph.create_node(f"{activity_id}_ShearVector", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, t_vector))
        self.graph.add((t_vector, CMSO.hasComponent_x, Literal(translation_vector[0], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_y, Literal(translation_vector[1], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_z, Literal(translation_vector[2], datatype=XSD.float),))

        #if plane is provided, add that as well
        if plane is not None:
            plane = self.graph.create_node(f"{activity_id}_Plane", CMSO.Plane)
            plane_vector = self.graph.create_node(f"{activity_id}_PlaneVector", CMSO.NormalVector)
            self.graph.add((activity, UNSAFECMSO.hasPlane, plane))
            self.graph.add((plane, CMSO.hasNormalVector, plane_vector))
            self.graph.add((plane_vector, CMSO.hasComponent_x, Literal(plane[0], datatype=XSD.float),))
            self.graph.add((plane_vector, CMSO.hasComponent_y, Literal(plane[1], datatype=XSD.float),))
            self.graph.add((plane_vector, CMSO.hasComponent_z, Literal(plane[2], datatype=XSD.float),))
            self.graph.add((plane, CMSO.hasDistanceFromOrigin, Literal(distance, datatype=XSD.float)))
        if original_sample != self.sample:
            sys.graph.add((sys.sample, PROV.wasDerivedFrom, original_sample))

    def copy_defects(self, parent_sample):
        if self.sample is None:
            return
        if parent_sample is None:
            return
        self.graph.copy_defects(self.sample, parent_sample)

    def plot3d(self, *args, **kwargs):
        try:
            from pyiron_atomistics.atomistics.structure.atoms import (
                ase_to_pyiron,
                pyiron_to_ase,
            )
        except ImportError:
            raise ImportError("Please install pyiron_atomistics")
        ase_structure = self.write.ase()
        pyiron_structure = ase_to_pyiron(ase_structure)
        return pyiron_structure.plot3d(*args, **kwargs)
