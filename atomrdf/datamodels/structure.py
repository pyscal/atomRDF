"""
Module defines the basic structure of atomic scale samples, including materials, crystal structures, unit cells, and simulation cells.
It also includes the definition of atom attributes and various types of defects.
"""

from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field, SkipValidation
from atomrdf.datamodels.basemodels import (
    TemplateMixin,
    RDFMixin,
    BaseModel,
)
from atomrdf.datamodels.workflow.property import Property as DataProperty
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
    DCAT,
)
import atomrdf.json_io as json_io
import atomrdf.datamodels.defects as defects
import atomrdf.datamodels.structure_io as structure_io
from atomrdf.utils import get_material, get_sample_id, get_sample_object, toPython
from atomrdf.datamodels.workflow.property import *
import atomrdf.properties as ap


# read element data file
file_location = os.path.dirname(__file__).split("/")
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__), "../data/element.yml")
with open(file_location, "r") as fin:
    element_identifiers = yaml.safe_load(fin)


class UnitCell(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    bravais_lattice: Optional[str] = None  # Ontology URI
    lattice_parameter: Optional[List[float]] = None  # Angstroms (implicit)
    angle: Optional[List[float]] = None  # Degrees (implicit)

    def to_graph(self, graph, sample_id, crystal_structure):
        unit_cell = graph.create_node(f"{sample_id}_UnitCell", CMSO.UnitCell)
        graph.add((crystal_structure, CMSO.hasUnitCell, unit_cell))

        if self.bravais_lattice is not None:
            bv = graph.create_node(
                f"{sample_id}_BravaisLattice", URIRef(self.bravais_lattice)
            )
            graph.add(
                (
                    unit_cell,
                    CMSO.hasBravaisLattice,
                    bv,
                )
            )

        if self.lattice_parameter is not None:
            lattice_parameter = graph.create_node(
                f"{sample_id}_LatticeParameter", CMSO.LatticeParameter
            )
            graph.add(
                (
                    unit_cell,
                    CMSO.hasLength_x,
                    Literal(self.lattice_parameter[0], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    unit_cell,
                    CMSO.hasLength_y,
                    Literal(self.lattice_parameter[1], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    unit_cell,
                    CMSO.hasLength_z,
                    Literal(self.lattice_parameter[2], datatype=XSD.float),
                )
            )

        if self.angle is not None:
            lattice_angle = graph.create_node(
                f"{sample_id}_LatticeAngle", CMSO.LatticeAngle
            )
            graph.add((unit_cell, CMSO.hasAngle, lattice_angle))
            graph.add(
                (
                    lattice_angle,
                    CMSO.hasAngle_alpha,
                    Literal(self.angle[0], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    lattice_angle,
                    CMSO.hasAngle_beta,
                    Literal(self.angle[1], datatype=XSD.float),
                )
            )
            graph.add(
                (
                    lattice_angle,
                    CMSO.hasAngle_gamma,
                    Literal(self.angle[2], datatype=XSD.float),
                )
            )

    @classmethod
    def from_graph(cls, graph, crystal_structure):
        unit_cell = graph.value(crystal_structure, CMSO.hasUnitCell)
        bv = graph.value(unit_cell, CMSO.hasBravaisLattice)
        bv = graph.value(bv, RDF.type)
        x = graph.value(unit_cell, CMSO.hasLength_x)
        y = graph.value(unit_cell, CMSO.hasLength_y)
        z = graph.value(unit_cell, CMSO.hasLength_z)
        angle = graph.value(unit_cell, CMSO.hasAngle)
        alpha = graph.value(angle, CMSO.hasAngle_alpha)
        beta = graph.value(angle, CMSO.hasAngle_beta)
        gamma = graph.value(angle, CMSO.hasAngle_gamma)
        datadict = {}
        if bv is not None:
            datadict["bravais_lattice"] = str(bv)
        if x is not None and y is not None and z is not None:
            datadict["lattice_parameter"] = [x.toPython(), y.toPython(), z.toPython()]
        if alpha is not None and beta is not None and gamma is not None:
            datadict["angle"] = [alpha.toPython(), beta.toPython(), gamma.toPython()]
        return cls(**datadict)


class CrystalStructure(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    name: Optional[str] = None
    spacegroup_symbol: Optional[str] = None
    spacegroup_number: Optional[int] = None
    unit_cell: Optional[UnitCell] = None

    def to_graph(self, graph, sample):
        # Custom logic to convert UnitCell to graph representation
        sample_id = get_sample_id(sample)
        material = get_material(graph, sample)

        crystal_structure = graph.create_node(
            f"{sample_id}_CrystalStructure", CMSO.CrystalStructure
        )
        graph.add((material, CMSO.hasStructure, crystal_structure))

        if self.spacegroup_symbol is not None:
            graph.add(
                (
                    crystal_structure,
                    CMSO.hasSpaceGroupSymbol,
                    Literal(self.spacegroup_symbol, datatype=XSD.string),
                )
            )
        if self.spacegroup_number is not None:
            graph.add(
                (
                    crystal_structure,
                    CMSO.hasSpaceGroupNumber,
                    Literal(self.spacegroup_number, datatype=XSD.integer),
                )
            )

        if self.unit_cell is not None:
            self.unit_cell.to_graph(graph, sample_id, crystal_structure)

    @classmethod
    def from_graph(cls, graph, sample):
        sample_id = get_sample_id(sample)
        material = get_material(graph, sample)
        crystal_structure = graph.value(material, CMSO.hasStructure)
        spacegroup_symbol = graph.value(crystal_structure, CMSO.hasSpaceGroupSymbol)
        spacegroup_number = graph.value(crystal_structure, CMSO.hasSpaceGroupNumber)
        unit_cell = UnitCell.from_graph(graph, crystal_structure)

        datadict = {
            "spacegroup_symbol": spacegroup_symbol,
            "spacegroup_number": spacegroup_number,
            "unit_cell": unit_cell,
        }
        return cls(**datadict)


class Material(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    element_ratio: Optional[dict] = None
    crystal_structure: Optional[CrystalStructure] = None

    def to_graph(self, graph, sample):
        sample_id = get_sample_id(sample)
        material = graph.create_node(f"{sample_id}_Material", CMSO.CrystallineMaterial)
        graph.add((sample, CMSO.hasMaterial, material))

        composition = self.element_ratio
        valid = True
        for e, r in composition.items():
            if e not in element_identifiers.keys():
                valid = False
                break

        if valid:
            chemical_species = graph.create_node(
                f"{sample_id}_ChemicalSpecies", CMSO.ChemicalSpecies
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
                    graph.add(
                        (element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float))
                    )

        if self.crystal_structure is not None:
            self.crystal_structure.to_graph(graph, sample)

    @classmethod
    def from_graph(cls, graph, sample):
        material = get_material(graph, sample)
        element_ratio = {}
        chemical_species = graph.value(sample, CMSO.hasSpecies)

        for element in graph.objects(chemical_species, CMSO.hasElement):
            symbol = graph.value(element, CMSO.hasChemicalSymbol)
            ratio = graph.value(element, CMSO.hasElementRatio)
            element_ratio[str(symbol)] = float(ratio)

        crystal_structure = CrystalStructure.from_graph(graph, sample)

        datadict = {
            "element_ratio": element_ratio,
            "crystal_structure": crystal_structure,
        }
        return cls(**datadict)


class SimulationCell(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    volume: Optional[DataProperty[float]] = None  # Keep - has unit (Angstrom³)
    number_of_atoms: Optional[int] = None
    length: Optional[List[float]] = None  # Angstroms (implicit)
    vector: Optional[List[List[float]]] = None  # Angstroms (implicit)
    angle: Optional[List[float]] = None  # Degrees (implicit)
    repetitions: Optional[List[int]] = None
    grain_size: Optional[float] = None  # angstroms (implicit)
    number_of_grains: Optional[int] = None

    def to_graph(self, graph, sample):
        sample_id = get_sample_id(sample)
        simulation_cell = graph.create_node(
            f"{sample_id}_SimulationCell", CMSO.SimulationCell
        )
        graph.add((sample, CMSO.hasSimulationCell, simulation_cell))

        # Only add volume if it exists and has a value
        if self.volume is not None and self.volume.value is not None:
            volume = graph.create_node(
                f"{sample_id}_Volume", ASMO.Volume, label="SimulationCellVolume"
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

        # Only add number of atoms if it's not None
        if self.number_of_atoms is not None:
            graph.add(
                (
                    sample,
                    CMSO.hasNumberOfAtoms,
                    Literal(self.number_of_atoms, datatype=XSD.integer),
                )
            )

        repetitions = self.repetitions
        if repetitions is not None:
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
        # Only add length if it exists and has data
        if self.length is not None and len(self.length) >= 3:
            simulation_cell_length = graph.create_node(
                f"{sample_id}_SimulationCellLength", CMSO.SimulationCellLength
            )
            graph.add((simulation_cell, CMSO.hasLength, simulation_cell_length))
            data = self.length
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

        # Only add vectors if they exist and have data
        if self.vector is not None and len(self.vector) >= 3:
            simulation_cell_vector_01 = graph.create_node(
                f"{sample_id}_SimulationCellVector_1", CMSO.SimulationCellVector
            )
            data = self.vector
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
                f"{sample_id}_SimulationCellVector_2", CMSO.SimulationCellVector
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
                f"{sample_id}_SimulationCellVector_3", CMSO.SimulationCellVector
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

        # Only add angles if they exist and have data
        if self.angle is not None and len(self.angle) >= 3:
            simulation_cell_angle = graph.create_node(
                f"{sample_id}_SimulationCellAngle", CMSO.SimulationCellAngle
            )
            data = self.angle
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

        if self.grain_size is not None:
            graph.add(
                (
                    simulation_cell,
                    CMSO.hasGrainSize,
                    Literal(self.grain_size, datatype=XSD.float),
                )
            )
        if self.number_of_grains is not None:
            graph.add(
                (
                    simulation_cell,
                    CMSO.hasNumberOfGrains,
                    Literal(self.number_of_grains, datatype=XSD.integer),
                )
            )

    @classmethod
    def from_graph(cls, graph, sample):
        simulation_cell = graph.value(sample, CMSO.hasSimulationCell)
        volume_item = graph.value(simulation_cell, CMSO.hasVolume)
        volume = graph.value(volume_item, ASMO.hasValue)
        number_of_atoms = graph.value(sample, CMSO.hasNumberOfAtoms)

        rx = graph.value(simulation_cell, CMSO.hasRepetition_x) or 1
        ry = graph.value(simulation_cell, CMSO.hasRepetition_y) or 1
        rz = graph.value(simulation_cell, CMSO.hasRepetition_z) or 1

        repetitions = [
            int(rx),
            int(ry),
            int(rz),
        ]
        simulation_cell_length = graph.value(
            simulation_cell,
            CMSO.hasLength,
        )
        length = None
        if simulation_cell_length is not None:
            lx = graph.value(simulation_cell_length, CMSO.hasLength_x)
            ly = graph.value(simulation_cell_length, CMSO.hasLength_y)
            lz = graph.value(simulation_cell_length, CMSO.hasLength_z)
            if lx is not None and ly is not None and lz is not None:
                length = [lx.toPython(), ly.toPython(), lz.toPython()]

        vector = []
        for v in graph.objects(simulation_cell, CMSO.hasVector):
            vector.append(
                [
                    toPython(graph.value(v, CMSO.hasComponent_x)),
                    toPython(graph.value(v, CMSO.hasComponent_y)),
                    toPython(graph.value(v, CMSO.hasComponent_z)),
                ]
            )
        cell_angle = graph.value(simulation_cell, CMSO.hasAngle)
        angle = [
            toPython(graph.value(cell_angle, CMSO.hasAngle_alpha)),
            toPython(graph.value(cell_angle, CMSO.hasAngle_beta)),
            toPython(graph.value(cell_angle, CMSO.hasAngle_gamma)),
        ]

        datadict = {
            "volume": {
                "value": volume,
            },
            "number_of_atoms": (
                int(number_of_atoms) if number_of_atoms is not None else None
            ),
            "repetitions": repetitions,
            "length": length,
            "vector": vector,
            "angle": angle,
        }
        return cls(**datadict)


class AtomAttribute(BaseModel, TemplateMixin):
    pid: Optional[str] = str(CMSO.AtomAttribute)
    position: SkipValidation[Optional[List[List[float]]]] = (
        None  # Angstroms (implicit) - validation skipped for performance
    )
    species: SkipValidation[Optional[List[str]]] = (
        None  # Chemical symbols - validation skipped for performance
    )

    def write_attributes(
        self, graph, sample_id, position_identifier, species_identifier
    ):
        datadict = {
            position_identifier: {
                "value": self.position,
                "label": "position",
            },
            species_identifier: {
                "value": self.species,
                "label": "species",
            },
        }
        outfile = os.path.join(graph.structure_store, str(sample_id).split(":")[-1])
        json_io.write_file(outfile, datadict)
        return os.path.relpath(outfile + ".json")

    def to_graph(self, graph, sample):
        # now we write out file
        sample_id = get_sample_id(sample)
        position_identifier = str(uuid.uuid4())
        species_identifier = str(uuid.uuid4())

        outfile = self.write_attributes(
            graph, sample_id, position_identifier, species_identifier
        )

        if self.position is not None:
            position = graph.create_node(f"{sample_id}_Position", CMSO.AtomAttribute)
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

        if self.species is not None:
            species = graph.create_node(f"{sample_id}_Species", CMSO.AtomAttribute)
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

    @classmethod
    def from_graph(cls, graph, sample):
        sample_id = get_sample_id(sample)

        # cell_vectors
        filepath = graph.value(URIRef(f"{sample_id}_Position"), CMSO.hasPath).toPython()
        position_identifier = graph.value(
            URIRef(f"{sample_id}_Position"), CMSO.hasIdentifier
        ).toPython()
        species_identifier = graph.value(
            URIRef(f"{sample_id}_Species"), CMSO.hasIdentifier
        ).toPython()

        # Resolve via structure_store so paths stored on a different machine/CWD still work
        store_dir = getattr(graph, "structure_store", None)
        if store_dir is not None:
            filepath = os.path.join(store_dir, os.path.basename(filepath))

        # open the file for reading
        with open(filepath, "r") as fin:
            data = json.load(fin)
            positions = data[position_identifier]["value"]
            species = data[species_identifier]["value"]

        return cls(
            position=positions,
            species=species,
        )


class AtomicScaleSample(BaseModel, TemplateMixin):
    pid: Optional[str] = str(CMSO.AtomicScaleSample)
    material: Optional[Material] = None
    simulation_cell: Optional[SimulationCell] = None
    atom_attribute: Optional[AtomAttribute] = None

    # add defects, all optional of course
    # point defects
    point_defect: Optional[defects.PointDefect] = None
    vacancy: Optional[defects.Vacancy] = None
    substitutional: Optional[defects.Substitutional] = None
    interstitial: Optional[defects.Interstitial] = None

    # dislocations
    dislocation: Optional[defects.Dislocation] = None
    edge_dislocation: Optional[defects.EdgeDislocation] = None
    screw_dislocation: Optional[defects.ScrewDislocation] = None
    mixed_dislocation: Optional[defects.MixedDislocation] = None

    # stacking faults
    stacking_fault: Optional[defects.StackingFault] = None

    # grain boundaries
    grain_boundary: Optional[defects.GrainBoundary] = None
    tilt_grain_boundary: Optional[defects.TiltGrainBoundary] = None
    twist_grain_boundary: Optional[defects.TwistGrainBoundary] = None
    symmetric_tilt_grain_boundary: Optional[defects.SymmetricalTiltGrainBoundary] = None
    mixed_grain_boundary: Optional[defects.MixedGrainBoundary] = None

    # properties
    calculated_property: Optional[List[CalculatedProperty]] = Field(
        default=[], description="Calculated properties from the simulation"
    )

    # defect complex
    defect_complex: Optional[defects.DefectComplex] = None

    def __init__(self, **data):
        super().__init__(**data)
        # Initialize private accumulator (not a Pydantic field)
        self._defect_complex_ids = []

    def to_graph_calculated_properties(self, graph):
        if self.calculated_property:
            for param in self.calculated_property:
                param_uri = param.to_graph(graph)
                graph.add(
                    (
                        URIRef(self.id),
                        ASMO.hasCalculatedProperty,
                        param_uri,
                    )
                )

    def from_graph_calculated_properties(cls, graph, sample_id):
        properties = []
        for prop_uri in graph.objects(URIRef(sample_id), ASMO.hasCalculatedProperty):
            prop = CalculatedProperty.from_graph(graph, prop_uri)
            properties.append(prop)
        cls.calculated_property = properties
        return cls

    def to_graph(self, graph, force=False):
        # if force - creates a new ID and saves the structure again
        if not force and self.id is not None:
            return self.id

        # the rest of the function is only if id isnt there or force is true
        name = f"sample:{str(uuid.uuid4())}"
        self.id = name
        sample = graph.create_node(name, CMSO.AtomicScaleSample, label=self.label)
        if self.material is not None:
            self.material.to_graph(graph, sample)
        if self.simulation_cell is not None:
            self.simulation_cell.to_graph(graph, sample)
        if self.atom_attribute is not None:
            self.atom_attribute.to_graph(
                graph,
                sample,
            )

        # now add calculated properties
        self.to_graph_calculated_properties(graph)

        # now call defect methods
        # Defects
        defect_fields = [
            "point_defect",
            "vacancy",
            "substitutional",
            "interstitial",
            "dislocation",
            "edge_dislocation",
            "screw_dislocation",
            "mixed_dislocation",
            "stacking_fault",
            "grain_boundary",
            "tilt_grain_boundary",
            "twist_grain_boundary",
            "symmetric_tilt_grain_boundary",
            "mixed_grain_boundary",
        ]

        # Reset the defect complex IDs accumulator
        self._defect_complex_ids = []

        # Map defect field names to their class names for URI construction
        defect_class_names = {
            "point_defect": "PointDefect",
            "vacancy": "Vacancy",
            "substitutional": "Substitutional",
            "interstitial": "Interstitial",
            "dislocation": "Dislocation",
            "edge_dislocation": "EdgeDislocation",
            "screw_dislocation": "ScrewDislocation",
            "mixed_dislocation": "MixedDislocation",
            "stacking_fault": "StackingFault",
            "grain_boundary": "GrainBoundary",
            "tilt_grain_boundary": "TiltGrainBoundary",
            "twist_grain_boundary": "TwistGrainBoundary",
            "symmetric_tilt_grain_boundary": "SymmetricTiltGrainBoundary",
            "mixed_grain_boundary": "MixedGrainBoundary",
        }

        for defect in defect_fields:
            obj = getattr(self, defect, None)
            if obj is not None:
                if isinstance(obj, BaseModel) and obj.model_fields_set:
                    if hasattr(obj, "to_graph"):
                        obj.to_graph(graph, sample)

                        # If this defect belongs to defect_complex, construct and collect its URI
                        if self.defect_complex and defect in self.defect_complex.ids:
                            # Construct the defect URI (matches the pattern used in defect.to_graph())
                            defect_uri = f"{sample}_{defect_class_names[defect]}"
                            self._defect_complex_ids.append(defect_uri)

        # Now serialize DefectComplex with collected defect URIs
        if self.defect_complex is not None:
            self.defect_complex.to_graph(graph, sample, self._defect_complex_ids)

        # Add content hash to the graph for deduplication (skip validation for external vocab)
        content_hash = self._compute_hash()
        graph.add((sample, DCAT.checksum, Literal(content_hash, datatype=XSD.string)))

        return self.id

    @classmethod
    def from_graph(cls, graph, sample_id):
        from typing import get_origin, get_args

        kwargs = {}
        sample = get_sample_object(sample_id)

        # try a type query first
        sample_type = graph.value(sample, RDF.type)
        if sample_type is None:
            raise ValueError(f"Sample {sample_id} not found in graph.")

        # material, simulation_cell, atom_attribute handled separately (if needed)
        kwargs["material"] = Material.from_graph(graph, sample)
        kwargs["simulation_cell"] = SimulationCell.from_graph(graph, sample)
        kwargs["atom_attribute"] = AtomAttribute.from_graph(graph, sample)

        defect_fields = [
            "point_defect",
            "vacancy",
            "substitutional",
            "interstitial",
            "dislocation",
            "edge_dislocation",
            "screw_dislocation",
            "mixed_dislocation",
            "stacking_fault",
            "grain_boundary",
            "tilt_grain_boundary",
            "twist_grain_boundary",
            "symmetric_tilt_grain_boundary",
            "mixed_grain_boundary",
        ]

        # Loop over defect fields
        for field in defect_fields:
            field_type = cls.model_fields[field].annotation

            # Unwrap Optional types (Union[X, None] -> X)
            origin = get_origin(field_type)
            if origin is not None:
                # Check if it's Optional/Union
                args = get_args(field_type)
                if len(args) == 2 and type(None) in args:
                    # It's Optional[X], extract X
                    field_type = args[0] if args[0] is not type(None) else args[1]

            if hasattr(field_type, "from_graph"):
                try:
                    result = field_type.from_graph(graph, sample)
                    kwargs[field] = result
                except Exception as e:
                    kwargs[field] = None

        cls = cls(**kwargs)
        cls.id = sample_id
        cls = cls.from_graph_calculated_properties(graph, sample_id)
        return cls

    @classmethod
    def from_file(
        cls,
        filename,
        format="lammps-dump",
        species=None,
        lattice=None,
        lattice_constant=None,
        basis_box=None,
        basis_positions=None,
        repeat=None,
        graph=None,
    ):
        """
        Read structure from file and create an AtomicScaleSample instance.

        Parameters
        ----------
        filename : str
            Path to the structure file
        format : str, optional
            File format (default: 'lammps-dump'). Any format supported by ASE.
        species : list, optional
            If provided, LAMMPS types will be matched to species. For example, if types 1 and 2 exist
            in the input file, and species = ['Li', 'Al'] is given, type 1 will be matched to 'Li' and
            type 2 will be matched to 'Al'
        lattice : str, optional
            Crystal structure name (e.g., 'bcc', 'fcc', 'hcp', 'diamond', 'l12', 'b2').
            If provided, metadata such as unit cell, space group, etc. are automatically added.
        lattice_constant : float, optional
            Lattice constant of the system
        basis_box : list of lists, optional
            3x3 matrix specifying the basis unit cell. Not required if lattice is provided.
        basis_positions : list of lists, optional
            Nx3 array specifying relative positions of atoms in the unit cell.
            Not required if lattice is provided.
        repeat : tuple or int, optional
            Number of repetitions of the unit cell in each direction.
        graph : KnowledgeGraph, optional
            If provided, the structure will be added to the graph.

        Returns
        -------
        AtomicScaleSample
            The created sample instance

        Examples
        --------
        >>> sample = AtomicScaleSample.from_file('structure.lmp', format='lammps-dump')
        >>> sample = AtomicScaleSample.from_file('POSCAR', format='vasp',
        ...                                       lattice='bcc', lattice_constant=2.87)
        """
        from ase.io import read as ase_read
        from pyscal3.core import structure_dict
        from atomrdf.build.bulk import _generate_atomic_sample_data
        from atomrdf.build.buildutils import _declass

        # Read structure with ASE
        atoms = ase_read(filename, format=format)

        # Handle species mapping for LAMMPS
        if species is not None:
            types = atoms.get_array("type") if "type" in atoms.arrays else None
            if types is not None:
                new_symbols = [species[int(t) - 1] for t in types]
                atoms.set_chemical_symbols(new_symbols)

        # Build metadata dict (sdict) for supplementing structure information
        sdict = {}

        if lattice is not None:
            # Get structure information from known lattices
            from pyscal3.core import structure_dict

            if lattice in structure_dict.keys():
                lattice_info = structure_dict[lattice].get("conventional", {})
                sdict["structure"] = lattice

                # Try to get spacegroup info
                if "spacegroup_symbol" in lattice_info:
                    sdict["spacegroup_symbol"] = lattice_info["spacegroup_symbol"]
                if "spacegroup_number" in lattice_info:
                    sdict["spacegroup_number"] = lattice_info["spacegroup_number"]

        if lattice_constant is not None:
            lattice_constant = _declass(lattice_constant)
            sdict["a"] = lattice_constant
            sdict["b"] = lattice_constant
            sdict["c"] = lattice_constant

        # If no lattice info provided, try to determine from structure
        if not sdict and lattice is None:
            try:
                # Try to get spacegroup info from the atoms object
                spacegroup_symbol = ap.get_spacegroup_symbol(atoms)
                spacegroup_number = ap.get_spacegroup_number(atoms)
                if spacegroup_symbol:
                    sdict["spacegroup_symbol"] = spacegroup_symbol
                if spacegroup_number:
                    sdict["spacegroup_number"] = spacegroup_number
            except:
                # If we can't determine spacegroup, that's okay
                pass

        # Generate sample data using the helper function from bulk.py
        data = _generate_atomic_sample_data(atoms, sdict if sdict else None, repeat)

        # Create the AtomicScaleSample instance
        sample = cls(**data)

        # Optionally add to graph
        if graph is not None:
            sample.to_graph(graph)
            atoms.info["id"] = sample.id
            atoms.info["graph"] = graph

        return sample

    @classmethod
    def from_repository(
        cls,
        repository="materials_project",
        api_key=None,
        material_ids=None,
        chemical_system=None,
        is_stable=True,
        conventional=True,
        graph=None,
    ):
        """
        Fetch structure(s) from an external repository and create AtomicScaleSample instance(s).

        Parameters
        ----------
        repository : str, optional
            Repository name. Currently supports: 'materials_project' (default).
        api_key : str
            API key for the repository.
        material_ids : list of str, optional
            List of material IDs to fetch. For Materials Project, these are mp-ids like ['mp-149', 'mp-13'].
        chemical_system : str, optional
            Chemical system string (e.g., 'Fe-C', 'Li-Co-O'). If provided, all stable materials
            in this system will be fetched.
        is_stable : bool, optional
            If True (default), only fetch stable materials. Only used with chemical_system.
        conventional : bool, optional
            If True (default), use conventional cell. If False, use primitive cell.
        graph : KnowledgeGraph, optional
            If provided, the structure(s) will be added to the graph.

        Returns
        -------
        AtomicScaleSample or list of AtomicScaleSample
            If a single material is fetched, returns AtomicScaleSample.
            If multiple materials are fetched, returns a list of AtomicScaleSample instances.

        Raises
        ------
        ValueError
            If neither material_ids nor chemical_system is provided.
        ImportError
            If the required repository client library is not installed.

        Examples
        --------
        Fetch a single material by ID:
        >>> sample = AtomicScaleSample.from_repository(
        ...     repository='materials_project',
        ...     api_key='your_api_key',
        ...     material_ids=['mp-149']
        ... )

        Fetch all stable materials in a chemical system:
        >>> samples = AtomicScaleSample.from_repository(
        ...     repository='materials_project',
        ...     api_key='your_api_key',
        ...     chemical_system='Fe-C'
        ... )

        Fetch and add to graph:
        >>> kg = KnowledgeGraph()
        >>> sample = AtomicScaleSample.from_repository(
        ...     repository='materials_project',
        ...     api_key='your_api_key',
        ...     material_ids=['mp-149'],
        ...     graph=kg
        ... )
        """
        if repository.lower() == "materials_project":
            return cls._from_materials_project(
                api_key=api_key,
                material_ids=material_ids,
                chemical_system=chemical_system,
                is_stable=is_stable,
                conventional=conventional,
                graph=graph,
            )
        else:
            raise ValueError(
                f"Repository '{repository}' is not supported. "
                "Currently supported: 'materials_project'"
            )

    @classmethod
    def _from_materials_project(
        cls,
        api_key,
        material_ids=None,
        chemical_system=None,
        is_stable=True,
        conventional=True,
        graph=None,
    ):
        """
        Internal method to fetch structures from Materials Project.
        """
        try:
            from mp_api.client import MPRester
        except ImportError:
            raise ImportError(
                "mp-api is not installed. Please install it with: pip install mp-api"
            )

        if api_key is None:
            raise ValueError("api_key is required for Materials Project")

        if (chemical_system is None) and (material_ids is None):
            raise ValueError(
                "Please provide either 'chemical_system' or 'material_ids'"
            )

        rest = {
            "use_document_model": False,
            "include_user_agent": True,
            "api_key": api_key,
        }

        with MPRester(**rest) as mpr:
            if chemical_system is not None:
                docs = mpr.materials.summary.search(
                    chemsys=chemical_system, is_stable=is_stable
                )
            else:
                docs = mpr.materials.summary.search(material_ids=material_ids)

        # Process documents and create samples
        samples = []

        for doc in docs:
            struct = doc["structure"]
            if conventional:
                aseatoms = struct.to_conventional().to_ase_atoms()
            else:
                aseatoms = struct.to_primitive().to_ase_atoms()

            symmetry = doc["symmetry"]

            # Generate sample data
            from atomrdf.build.bulk import _generate_atomic_sample_data

            data = _generate_atomic_sample_data(aseatoms)
            sample = cls(**data)

            # Update spacegroup information
            sample.material.crystal_structure.spacegroup_symbol = symmetry["symbol"]
            sample.material.crystal_structure.spacegroup_number = symmetry["number"]

            # Add energy as a calculated property
            if "energy_per_atom" in doc and doc["energy_per_atom"] is not None:
                from atomrdf.datamodels.workflow.property import CalculatedProperty

                energy_prop = CalculatedProperty(
                    basename="PotentialEnergy",
                    label="Potential energy per atom",
                    value=float(doc["energy_per_atom"]),
                    unit="EV",
                )
                sample.calculated_property = [energy_prop]

            # Add to graph if provided
            if graph is not None:
                sample.to_graph(graph)
                aseatoms.info["id"] = sample.id
                aseatoms.info["graph"] = graph

            samples.append(sample)

        # Return single sample or list
        if len(samples) == 1:
            return samples[0]
        else:
            return samples

    def update_attributes(self, atoms, repeat=None):
        """
        Update the atom attributes based on the provided ASE Atoms object.
        This would also reset the id, since the structure has changed.
        """
        self.id = None
        self.material.element_ratio = ap.get_chemical_composition(atoms)
        self.simulation_cell.volume.value = ap.get_cell_volume(atoms)
        self.simulation_cell.number_of_atoms = ap.get_number_of_atoms(atoms)
        self.simulation_cell.length = ap.get_simulation_cell_length(atoms)
        self.simulation_cell.vector = ap.get_simulation_cell_vector(atoms)
        self.simulation_cell.angle = ap.get_simulation_cell_angle(atoms)
        if repeat is not None:
            if isinstance(repeat, int):
                self.simulation_cell.repetitions = (repeat, repeat, repeat)
            else:
                self.simulation_cell.repetitions = repeat

        self.atom_attribute.position = atoms.get_positions().tolist()
        self.atom_attribute.species = atoms.get_chemical_symbols()

    def to_structure(self, format="ase"):
        if format == "ase":
            return structure_io.sample_to_ase(self)
        else:
            raise ValueError(f"Unsupported format: {format}")

    def to_file(
        self,
        outfile,
        format,
        copy_from=None,
        pseudo_files=None,
    ):
        """
        Write the structure to a file in the specified format.

        Parameters
        ----------
        outfile : str
            The path to the output file.
        format : str, optional
            The format of the output file. Defaults to 'lammps-dump'.
        copy_from : str, optional
            If provided, input options for quantum-espresso format will be copied from
            the given file. Structure specific information will be replaced.
            Note that the validity of input file is not checked.
        pseudo_files : list, optional
            if provided, add the pseudopotential filenames to file.
            Should be in alphabetical order of chemical species symbols.

        Returns
        -------
        None
        """
        structure_io.write(
            self,
            outfile,
            format=format,
            copy_from=copy_from,
            pseudo_files=pseudo_files,
        )

    def _compute_hash(self, precision=6):
        """
        Compute a deterministic hash of the sample structure (internal method).

        Excludes 'id' and 'graph' fields and rounds floating-point values
        to the specified precision to ensure consistent hashing across
        equivalent structures.

        Parameters
        ----------
        precision : int, optional
            Number of decimal places for rounding floats (default: 6)

        Returns
        -------
        str
            MD5 hash (hexadecimal string) of the sample content

        Examples
        --------
        >>> sample1 = AtomicScaleSample(**sample_data)
        >>> sample2 = AtomicScaleSample(**sample_data)
        >>> sample1._compute_hash() == sample2._compute_hash()
        True
        """
        import hashlib
        import json

        # Convert to dict and exclude id/label fields
        data = self.model_dump(exclude={"id", "label"})

        # Helper to round all floats in nested structures
        def round_floats(obj, decimals):
            if isinstance(obj, float):
                return round(obj, decimals)
            elif isinstance(obj, dict):
                return {k: round_floats(v, decimals) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [round_floats(v, decimals) for v in obj]
            else:
                return obj

        # Round all floats to avoid precision issues
        data_rounded = round_floats(data, precision)

        # Create deterministic JSON string (sorted keys)
        json_str = json.dumps(data_rounded, sort_keys=True)

        # Compute and return MD5 hash
        return hashlib.md5(json_str.encode()).hexdigest()
