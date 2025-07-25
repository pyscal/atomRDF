from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
from pydantic import BaseModel, Field
from atomrdf.datamodels.basemodels import TemplateMixin, DataProperty
from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.namespace import (
    CMSO,
    LDO,
    PLDO,
    PODO,
    CDCO,
    PROV,
    Literal,
    ASMO,
)
import atomrdf.json_io as json_io

# read element data file
file_location = os.path.dirname(__file__).split("/")
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__), "../data/element.yml")
with open(file_location, "r") as fin:
    element_identifiers = yaml.safe_load(fin)


class UnitCell(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    bravais_lattice: Optional[DataProperty[str]] = None
    lattice_parameter: Optional[DataProperty[List[float]]] = None
    angle: Optional[DataProperty[List[float]]] = None

    def to_graph(self, graph, name, crystal_structure):
        unit_cell = graph.create_node(f"{name}_UnitCell", CMSO.UnitCell)
        graph.add((crystal_structure, CMSO.hasUnitCell, unit_cell))

        bv = graph.create_node(
            f"{name}_BravaisLattice", URIRef(self.bravais_lattice.value)
        )
        graph.add(
            (
                unit_cell,
                CMSO.hasBravaisLattice,
                bv,
            )
        )
        lattice_parameter = graph.create_node(
            f"{name}_LatticeParameter", CMSO.LatticeParameter
        )
        graph.add(
            (
                unit_cell,
                CMSO.hasLength_x,
                Literal(self.lattice_parameter.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                unit_cell,
                CMSO.hasLength_y,
                Literal(self.lattice_parameter.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                unit_cell,
                CMSO.hasLength_z,
                Literal(self.lattice_parameter.value[2], datatype=XSD.float),
            )
        )

        lattice_angle = graph.create_node(f"{name}_LatticeAngle", CMSO.LatticeAngle)
        graph.add((unit_cell, CMSO.hasAngle, lattice_angle))
        graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_alpha,
                Literal(self.angle.value[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_beta,
                Literal(self.angle.value[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_gamma,
                Literal(self.angle.value[2], datatype=XSD.float),
            )
        )


class CrystalStructure(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    name: Optional[DataProperty[str]] = None
    spacegroup_symbol: Optional[DataProperty[str]] = None
    spacegroup_number: Optional[DataProperty[int]] = None
    unit_cell: Optional[UnitCell] = None

    def to_graph(self, graph, name):
        # Custom logic to convert UnitCell to graph representation
        crystal_structure = graph.create_node(
            f"{name}_CrystalStructure", CMSO.CrystalStructure
        )
        graph.add(
            (
                crystal_structure,
                CMSO.hasSpaceGroupSymbol,
                Literal(self.spacegroup_symbol.value, datatype=XSD.string),
            )
        )
        graph.add(
            (
                crystal_structure,
                CMSO.hasSpaceGroupNumber,
                Literal(self.spacegroup_number.value, datatype=XSD.integer),
            )
        )
        self.unit_cell.to_graph(graph, name, crystal_structure)


class Material(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    element_ratio: Optional[DataProperty[dict]] = None
    crystal_structure: Optional[CrystalStructure] = None

    def to_graph(self, graph, name, sample):
        material = graph.create_node(f"{name}_Material", CMSO.CrystallineMaterial)
        graph.add((sample, CMSO.hasMaterial, material))

        composition = self.element_ratio.value
        valid = True
        for e, r in composition.items():
            if e in element_identifiers.keys():
                valid = False
                break

        if valid:
            chemical_species = graph.create_node(
                f"{name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            graph.add((sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_identifiers.keys():
                    element = graph.create_node(
                        element_identifiers[e], CMSO.ChemicalElement
                    )
                    graph.add((chemical_species, CMSO.hasElement, element))
                    graph.add(
                        (
                            element,
                            CMSO.hasChemicalSymbol,
                            Literal(e, datatype=XSD.string),
                        )
                    )
                    self.graph.add(
                        (element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float))
                    )

        self.crystal_structure.to_graph(graph, name)


class SimulationCell(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    volume: Optional[DataProperty[float]] = None
    number_of_atoms: Optional[DataProperty[int]] = None
    length: Optional[DataProperty[List[float]]] = None
    vector: Optional[DataProperty[List[List[float]]]] = None
    angle: Optional[DataProperty[List[float]]] = None
    repetitions: Optional[DataProperty[List[int]]] = None

    def to_graph(self, graph, name, sample):
        simulation_cell = graph.create_node(
            f"{name}_SimulationCell", CMSO.SimulationCell
        )
        graph.add((sample, CMSO.hasSimulationCell, simulation_cell))
        volume = graph.create_node(
            f"{name}_Volume", ASMO.Volume, label="SimulationCellVolume"
        )
        graph.add((simulation_cell, CMSO.hasVolume, volume))
        graph.add(
            (
                volume,
                ASMO.hasValue,
                Literal(
                    np.round(self.volume.value, decimals=2),
                    datatype=XSD.float,
                ),
            )
        )
        graph.add(
            (
                volume,
                ASMO.hasUnit,
                URIRef(f"http://qudt.org/vocab/unit/ANGSTROM3"),
            )
        )
        graph.add(
            (
                sample,
                CMSO.hasNumberOfAtoms,
                Literal(self.number_of_atoms.value, datatype=XSD.integer),
            )
        )

        repetitions = self.repetitions.value
        graph.add(
            (
                simulation_cell,
                CMSO.hasRepetition_x,
                Literal(repetitions[0], datatype=XSD.integer),
            )
        )
        graph.add(
            (
                simulation_cell,
                CMSO.hasRepetition_y,
                Literal(repetitions[1], datatype=XSD.integer),
            )
        )
        graph.add(
            (
                simulation_cell,
                CMSO.hasRepetition_z,
                Literal(repetitions[2], datatype=XSD.integer),
            )
        )

        simulation_cell_length = graph.create_node(
            f"{name}_SimulationCellLength", CMSO.SimulationCellLength
        )
        graph.add((simulation_cell, CMSO.hasLength, simulation_cell_length))
        data = self.length.value
        graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_x,
                Literal(data[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_y,
                Literal(data[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_z,
                Literal(data[2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_01 = graph.create_node(
            f"{name}_SimulationCellVector_1", CMSO.SimulationCellVector
        )
        data = self.vector.value
        graph.add((simulation_cell, CMSO.hasVector, simulation_cell_vector_01))
        graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_x,
                Literal(data[0][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_y,
                Literal(data[0][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_z,
                Literal(data[0][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_02 = graph.create_node(
            f"{name}_SimulationCellVector_2", CMSO.SimulationCellVector
        )
        graph.add((simulation_cell, CMSO.hasVector, simulation_cell_vector_02))
        graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_x,
                Literal(data[1][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_y,
                Literal(data[1][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_z,
                Literal(data[1][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_03 = graph.create_node(
            f"{name}_SimulationCellVector_3", CMSO.SimulationCellVector
        )
        graph.add((simulation_cell, CMSO.hasVector, simulation_cell_vector_03))
        graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_x,
                Literal(data[2][0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_y,
                Literal(data[2][1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_z,
                Literal(data[2][2], datatype=XSD.float),
            )
        )

        simulation_cell_angle = graph.create_node(
            f"{name}_SimulationCellAngle", CMSO.SimulationCellAngle
        )
        data = self.angle.value
        graph.add((simulation_cell, CMSO.hasAngle, simulation_cell_angle))
        graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_alpha,
                Literal(data[0], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_beta,
                Literal(data[1], datatype=XSD.float),
            )
        )
        graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_gamma,
                Literal(data[2], datatype=XSD.float),
            )
        )


class AtomAttribute(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    position: Optional[DataProperty[List[List[float]]]] = None
    species: Optional[DataProperty[List[str]]] = None

    def write_attributes(self, graph, name, position_identifier, species_identifier):
        datadict = {
            position_identifier: {
                "value": self.position.value,
                "label": "position",
            },
            species_identifier: {
                "value": self.species.value,
                "label": "species",
            },
        }
        outfile = os.path.join(graph.structure_store, str(name).split(":")[-1])
        json_io.write_file(outfile, datadict)
        return os.path.relpath(outfile + ".json")

    def to_graph(self, graph, name, sample):
        # now we write out file
        position_identifier = str(uuid.uuid4())
        species_identifier = str(uuid.uuid4())

        outfile = self.write_attributes(
            graph, name, position_identifier, species_identifier
        )

        if self.position.value is not None:
            position = graph.create_node(f"{name}_Position", CMSO.AtomAttribute)
            graph.add(
                (
                    sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    position,
                )
            )
            graph.add(
                (position, CMSO.hasName, Literal("Position", datatype=XSD.string))
            )
            graph.add(
                (
                    position,
                    CMSO.hasIdentifier,
                    Literal(position_identifier, datatype=XSD.string),
                )
            )
            graph.add((position, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))

        if self.species.value is not None:
            species = graph.create_node(f"{name}_Species", CMSO.AtomAttribute)
            graph.add(
                (
                    sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    species,
                )
            )
            graph.add((species, CMSO.hasName, Literal("Species", datatype=XSD.string)))
            graph.add(
                (
                    species,
                    CMSO.hasIdentifier,
                    Literal(species_identifier, datatype=XSD.string),
                )
            )
            graph.add((species, CMSO.hasPath, Literal(outfile, datatype=XSD.string)))


class AtomicScaleSample(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    material: Optional[Material] = None
    simulation_cell: Optional[SimulationCell] = None
    atom_attribute: Optional[AtomAttribute] = None

    def to_graph(self, graph):
        name = f"sample:{str(uuid.uuid4())}"
        sample = graph.create_node(name, CMSO.AtomicScaleSample, label=self.label)
        self.material.to_graph(graph, name, sample)
        self.simulation_cell.to_graph(graph, name, sample)
        self.atom_attribute.to_graph(
            graph,
            name,
            sample,
        )
        self.id = name
