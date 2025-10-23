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
from pydantic import Field
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
)
import atomrdf.json_io as json_io
import atomrdf.datamodels.defects as defects
import atomrdf.datamodels.structure_io as structure_io
from atomrdf.utils import get_material, get_sample_id, get_sample_object, toPython
import atomrdf.properties as ap


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

    def to_graph(self, graph, sample_id, crystal_structure):
        unit_cell = graph.create_node(f"{sample_id}_UnitCell", CMSO.UnitCell)
        graph.add((crystal_structure, CMSO.hasUnitCell, unit_cell))

        if self.bravais_lattice.value is not None: 
            bv = graph.create_node(
                f"{sample_id}_BravaisLattice", URIRef(self.bravais_lattice.value)
            )
            graph.add(
                (
                    unit_cell,
                    CMSO.hasBravaisLattice,
                    bv,
                )
            )
        
        if self.lattice_parameter.value is not None:
            lattice_parameter = graph.create_node(
                f"{sample_id}_LatticeParameter", CMSO.LatticeParameter
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

        if self.angle.value is not None:
            lattice_angle = graph.create_node(
                f"{sample_id}_LatticeAngle", CMSO.LatticeAngle
            )
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
            datadict["bravais_lattice"] = {"value": str(bv)}
        if x is not None and y is not None and z is not None:
            datadict["lattice_parameter"] = {"value": [x.toPython(), 
                                                       y.toPython(), 
                                                       z.toPython()]}
        if alpha is not None and beta is not None and gamma is not None:
            datadict["angle"] = {"value": [alpha.toPython(), 
                                           beta.toPython(), 
                                           gamma.toPython()]}
        return cls(**datadict)


class CrystalStructure(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    name: Optional[DataProperty[str]] = None
    spacegroup_symbol: Optional[DataProperty[str]] = None
    spacegroup_number: Optional[DataProperty[int]] = None
    unit_cell: Optional[UnitCell] = None

    def to_graph(self, graph, sample):
        # Custom logic to convert UnitCell to graph representation
        sample_id = get_sample_id(sample)
        material = get_material(graph, sample)

        crystal_structure = graph.create_node(
            f"{sample_id}_CrystalStructure", CMSO.CrystalStructure
        )
        graph.add((material, CMSO.hasStructure, crystal_structure))

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
            "spacegroup_symbol": {"value": spacegroup_symbol},
            "spacegroup_number": {"value": spacegroup_number},
            "unit_cell": unit_cell,
        }
        return cls(**datadict)


class Material(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    element_ratio: Optional[DataProperty[dict]] = None
    crystal_structure: Optional[CrystalStructure] = None

    def to_graph(self, graph, sample):
        sample_id = get_sample_id(sample)
        material = graph.create_node(f"{sample_id}_Material", CMSO.CrystallineMaterial)
        graph.add((sample, CMSO.hasMaterial, material))

        composition = self.element_ratio.value
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
            "element_ratio": {"value": element_ratio},
            "crystal_structure": crystal_structure,
        }
        return cls(**datadict)


class SimulationCell(BaseModel, TemplateMixin):
    pid: Optional[str] = None
    volume: Optional[DataProperty[float]] = None
    number_of_atoms: Optional[DataProperty[int]] = None
    length: Optional[DataProperty[List[float]]] = None
    vector: Optional[DataProperty[List[List[float]]]] = None
    angle: Optional[DataProperty[List[float]]] = None
    repetitions: Optional[DataProperty[List[int]]]

    def to_graph(self, graph, sample):
        sample_id = get_sample_id(sample)
        simulation_cell = graph.create_node(
            f"{sample_id}_SimulationCell", CMSO.SimulationCell
        )
        graph.add((sample, CMSO.hasSimulationCell, simulation_cell))
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
        graph.add(
            (
                sample,
                CMSO.hasNumberOfAtoms,
                Literal(self.number_of_atoms.value, datatype=XSD.integer),
            )
        )

        repetitions = self.repetitions.value
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
        simulation_cell_length = graph.create_node(
            f"{sample_id}_SimulationCellLength", CMSO.SimulationCellLength
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
            f"{sample_id}_SimulationCellVector_1", CMSO.SimulationCellVector
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

        simulation_cell_angle = graph.create_node(
            f"{sample_id}_SimulationCellAngle", CMSO.SimulationCellAngle
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
            "volume": {"value": volume,},
            "number_of_atoms": {"value": int(number_of_atoms)},
            "repetitions": {"value": repetitions},
            "length": {"value": length},
            "vector": {"value": vector},
            "angle": {"value": angle},
        }
        return cls(**datadict)


class AtomAttribute(BaseModel, TemplateMixin):
    pid: Optional[str] = CMSO.AtomAttribute.uri
    position: Optional[DataProperty[List[List[float]]]] = None
    species: Optional[DataProperty[List[str]]] = None

    def write_attributes(
        self, graph, sample_id, position_identifier, species_identifier
    ):
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

        if self.position.value is not None:
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

        if self.species.value is not None:
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

        # open the file for reading
        with open(filepath, "r") as fin:
            data = json.load(fin)
            positions = data[position_identifier]["value"]
            species = data[species_identifier]["value"]
        
        position=DataProperty(value=positions)
        species=DataProperty(value=species)
        
        return cls(
            position=position,
            species=species,
        )


class AtomicScaleSample(BaseModel, TemplateMixin):
    pid: Optional[str] = CMSO.AtomicScaleSample.uri
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

    #properties


    def to_graph(self, graph, force=False):
        # if force - creates a new ID and saves the structure again
        if not force and self.id is not None:
            return

        # the rest of the function is only if id isnt there or force is true
        name = f"sample:{str(uuid.uuid4())}"
        self.id = name
        sample = graph.create_node(name, CMSO.AtomicScaleSample, label=self.label)
        self.material.to_graph(graph, sample)
        self.simulation_cell.to_graph(graph, sample)
        self.atom_attribute.to_graph(
            graph,
            sample,
        )

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

        for defect in defect_fields:
            obj = getattr(self, defect)
            if isinstance(obj, BaseModel) and obj.model_fields_set:
                if hasattr(obj, "to_graph"):
                    obj.to_graph(graph, sample)

    @classmethod
    def from_graph(cls, graph, sample_id):
        kwargs = {}
        sample = get_sample_object(sample_id)

        #try a type query first
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
            if hasattr(field_type, "from_graph"):
                try:
                    kwargs[field] = field_type.from_graph(graph, sample)
                except Exception:
                    kwargs[field] = None  # or skip/log

        cls = cls(**kwargs)
        cls.id = sample_id
        return cls

    def update_attributes(self, atoms, repeat=None):
        """
        Update the atom attributes based on the provided ASE Atoms object.
        This would also reset the id, since the structure has changed.
        """
        self.id = None
        self.material.element_ratio.value = ap.get_chemical_composition(atoms)
        self.simulation_cell.volume.value = ap.get_cell_volume(atoms)
        self.simulation_cell.number_of_atoms.value = ap.get_number_of_atoms(atoms)
        self.simulation_cell.length.value = ap.get_simulation_cell_length(atoms)
        self.simulation_cell.vector.value = ap.get_simulation_cell_vector(atoms)
        self.simulation_cell.angle.value = ap.get_simulation_cell_angle(atoms)
        if repeat is not None:
            if isinstance(repeat, int):
                self.simulation_cell.repetitions.value = (repeat, repeat, repeat)
            else:
                self.simulation_cell.repetitions.value = repeat

        self.atom_attribute.position.value = atoms.get_positions().tolist()
        self.atom_attribute.species.value = atoms.get_chemical_symbols()

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
