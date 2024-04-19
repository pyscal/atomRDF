"""
Workflows aspects for non-automated annotation of structures.

This consists of a workflow class which implements the necessary methods to serialise triples as needed.
Custom workflow solutions can be implemented. An example available here is pyiron.
The custom workflow env should implement the following functions:

_check_if_job_is_valid
_add_structure
_identify_method
extract_calculated_properties
inform_graph

See atomrdf.workflow.pyiron for more details
"""

from rdflib import Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef

import warnings
import numpy as np
import os
import copy
import ast
import uuid

from atomrdf.structure import System

# Move imports to another file
from atomrdf.namespace import PROV, CMSO, PODO, ASMO

# custom imports as needed
import atomrdf.workflow.pyiron as pi


class Workflow:
    def __init__(self, kg, environment="pyiron"):
        """
        Initialize the workflow environment

        Parameters
        ----------
        kg: pyscal-rdf KnowledgeGraph
        environment: string
            the workflow environment. This is used to import the necessary functions.

        """
        self.kg = kg
        if environment == "pyiron":
            self.wenv = pi
        else:
            raise ValueError("unknown workflow environment")

    def _prepare_job(self, workflow_object):

        self.wenv._check_if_job_is_valid(workflow_object)
        parent_structure, parent_sample, structure, sample = self.wenv._add_structures(
            workflow_object
        )
        method_dict = self.wenv._identify_method(workflow_object)

        if (structure is None) and (sample is None):
            raise ValueError("Either structure or sample should be specified")

        if sample is None:
            # its not added to graph yet
            structure.graph = self.kg
            structure.to_graph()
            sample = structure.sample

        if parent_sample is None:
            # its not added to graph yet
            if parent_structure is not None:
                parent_structure.graph = self.kg
                parent_structure.to_graph()
                parent_sample = parent_structure.sample

        self.structure = structure
        self.sample = sample
        self.mdict = method_dict
        self.main_id = method_dict["id"]
        self.parent_sample = parent_sample

    def _get_lattice_properties(
        self,
    ):
        if self.parent_sample is None:
            return

        material = list(
            [
                k[2]
                for k in self.kg.triples((self.parent_sample, CMSO.hasMaterial, None))
            ]
        )[0]
        crystal_structure = self.kg.value(material, CMSO.hasStructure)

        altname = self.kg.value(crystal_structure, CMSO.hasAltName)

        space_group = self.kg.value(crystal_structure, CMSO.hasSpaceGroup)
        space_group_symbol = self.kg.value(space_group, CMSO.hasSpaceGroupSymbol)
        space_group_number = self.kg.value(space_group, CMSO.hasSpaceGroupNumber)

        unit_cell = self.kg.value(crystal_structure, CMSO.hasUnitCell)
        blattice = self.kg.value(
            unit_cell,
            Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasBravaisLattice,
        )

        lattice_parameter = self.kg.value(unit_cell, CMSO.hasLatticeParameter)
        lattice_parameter_x = self.kg.value(lattice_parameter, CMSO.hasLength_x)
        lattice_parameter_y = self.kg.value(lattice_parameter, CMSO.hasLength_y)
        lattice_parameter_z = self.kg.value(lattice_parameter, CMSO.hasLength_z)

        lattice_angle = self.kg.value(unit_cell, CMSO.hasAngle)
        lattice_angle_x = self.kg.value(lattice_angle, CMSO.hasAngle_alpha)
        lattice_angle_y = self.kg.value(lattice_angle, CMSO.hasAngle_beta)
        lattice_angle_z = self.kg.value(lattice_angle, CMSO.hasAngle_gamma)

        targets = [
            altname,
            space_group_symbol,
            space_group_number,
            blattice,
            [lattice_parameter_x, lattice_parameter_y, lattice_parameter_z],
            [lattice_angle_x, lattice_angle_y, lattice_angle_z],
        ]

        self.structure._add_crystal_structure(targets=targets)

    def _add_inherited_properties(
        self,
    ):
        # Here we need to add inherited info: CalculatedProperties will be lost
        # Defects will be inherited
        if self.parent_sample is None:
            return

        parent_material = list(
            [
                k[2]
                for k in self.kg.triples((self.parent_sample, CMSO.hasMaterial, None))
            ]
        )[0]
        parent_defects = list(
            [x[2] for x in self.kg.triples((parent_material, CMSO.hasDefect, None))]
        )
        # now for each defect we copy add this to the final sample
        material = list(
            [k[2] for k in self.kg.triples((self.sample, CMSO.hasMaterial, None))]
        )[0]

        for defect in parent_defects:
            new_defect = URIRef(defect.toPython())
            self.kg.add((material, CMSO.hasDefect, new_defect))
            # now fetch all defect based info
            for triple in self.kg.triples((defect, None, None)):
                self.kg.add((new_defect, triple[1], triple[2]))

        # now add the special props for vacancy
        parent_simcell = self.kg.value(self.sample, CMSO.hasSimulationCell)
        simcell = self.kg.value(self.parent_sample, CMSO.hasSimulationCell)

        for triple in self.kg.triples(
            (parent_simcell, PODO.hasVacancyConcentration, None)
        ):
            self.kg.add((simcell, triple[1], triple[2]))
        for triple in self.kg.triples(
            (parent_simcell, PODO.hasNumberOfVacancies, None)
        ):
            self.kg.add((simcell, triple[1], triple[2]))

    def add_structural_relation(
        self,
    ):
        """
        Add structural relation between samples.

        This method adds the structural relation between the current sample and its parent sample.
        It also retrieves lattice properties and adds inherited properties, such as defect information.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.kg.add((self.sample, RDF.type, PROV.Entity))
        if self.parent_sample is not None:
            self.kg.add((self.parent_sample, RDF.type, PROV.Entity))
            self.kg.add((self.sample, PROV.wasDerivedFrom, self.parent_sample))
            self._get_lattice_properties()
            self._add_inherited_properties()

    def add_method(
            self,
        ):
        """
        Add the computational method and related information to the knowledge graph.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method adds the computational method and related information to the knowledge graph.
        It creates an activity node representing the method and adds it to the graph.
        The method is associated with the main activity using the `ASMO.hasComputationalMethod` property.
        The type of the method is determined based on the value of the `method` key in the `mdict` dictionary.
        The method-specific items are added to the graph based on the method type.
        The structure generation information is also added to the graph.

        """
        if self.mdict is None:
            return

        # add activity
        # ----------------------------------------------------------
        activity = URIRef(f"activity_{self.main_id}")
        self.kg.add((activity, RDF.type, PROV.Activity))

        # add method
        # ----------------------------------------------------------
        method = URIRef(f"method_{self.main_id}")
        if self.mdict["method"] == "MolecularStatics":
            self.kg.add((method, RDF.type, ASMO.MolecularStatics))
        elif self.mdict["method"] == "MolecularDynamics":
            self.kg.add((method, RDF.type, ASMO.MolecularDynamics))
        elif self.mdict["method"] == "DensityFunctionalTheory":
            self.kg.add((method, RDF.type, ASMO.DensityFunctionalTheory))
        self.kg.add((activity, ASMO.hasComputationalMethod, method))

        # choose if its rigid energy or structure optimisation
        # ----------------------------------------------------------
        if len(self.mdict["dof"]) == 0:
            self.kg.add(
                (
                    activity,
                    RDF.type,
                    Namespace(
                        "http://purls.helmholtz-metadaten.de/asmo/"
                    ).RigidEnergyCalculation,
                )
            )
        else:
            self.kg.add((activity, RDF.type, ASMO.StructureOptimization))
        # add DOFs
        for dof in self.mdict["dof"]:
            self.kg.add((activity, ASMO.hasRelaxationDOF, getattr(ASMO, dof)))

        # add method specific items
        if self.mdict["method"] in ["MolecularStatics", "MolecularDynamics"]:
            self._add_md(method, activity)
        elif self.mdict["method"] in ["DensityFunctionalTheory"]:
            self._add_dft(method, activity)

        # add that structure was generated
        self.kg.add((self.sample, PROV.wasGeneratedBy, activity))
        self._add_inputs(activity)
        self._add_outputs(activity)
        self._add_software(method)

    def to_graph(self, workflow_object):
        """
        Converts a workflow object to a graph representation.

        Parameters:
        - workflow_object: The workflow object to convert.

        Returns:
        - None
        """
        self._prepare_job(workflow_object)
        self.add_structural_relation()
        self.add_method()

    def _add_outputs(self, activity):
        if "outputs" in self.mdict.keys():
            for out in self.mdict["outputs"]:
                prop = self.kg.create_node(
                    f'{self.main_id}_{out["label"]}', CMSO.CalculatedProperty
                )
                self.kg.add((prop, RDFS.label, Literal(out["label"])))
                self.kg.add((prop, ASMO.hasValue, Literal(out["value"])))
                if "unit" in out.keys():
                    unit = out["unit"]
                    self.kg.add(
                        (
                            prop,
                            ASMO.hasUnit,
                            URIRef(f"http://qudt.org/vocab/unit/{unit}"),
                        )
                    )
                self.kg.add((prop, ASMO.wasCalculatedBy, activity))
                if out["associate_to_sample"]:
                    self.kg.add((self.sample, CMSO.hasCalculatedProperty, prop))

    def _add_inputs(self, activity):
        if "inputs" in self.mdict.keys():
            for inp in self.mdict["inputs"]:
                prop = self.kg.create_node(
                    f'{self.main_id}_{inp["label"]}', ASMO.InputParameter
                )
                self.kg.add((prop, RDFS.label, Literal(inp["label"])))
                self.kg.add((prop, ASMO.hasValue, Literal(inp["value"])))
                if "unit" in inp.keys():
                    unit = inp["unit"]
                    self.kg.add(
                        (
                            prop,
                            ASMO.hasUnit,
                            URIRef(f"http://qudt.org/vocab/unit/{unit}"),
                        )
                    )
                self.kg.add((activity, ASMO.hasInputParameter, prop))

    def _add_software(self, method):
        # finally add software
        wfagent = None
        if "workflow_manager" in self.mdict.keys():
            wfagent = self.kg.create_node(
                self.mdict["workflow_manager"]["uri"], PROV.SoftwareAgent
            )
            self.kg.add(
                (wfagent, RDFS.label, Literal(self.mdict["workflow_manager"]["label"]))
            )
            self.kg.add((method, PROV.wasAssociatedWith, wfagent))

        for software in self.mdict["software"]:
            agent = self.kg.create_node(software["uri"], PROV.SoftwareAgent)
            self.kg.add((agent, RDFS.label, Literal(software["label"])))
            if wfagent is not None:
                self.kg.add((wfagent, PROV.actedOnBehalfOf, agent))
            else:
                self.kg.add((method, PROV.wasAssociatedWith, agent))

    def _add_md(self, method, activity):
        self.kg.add(
            (method, ASMO.hasStatisticalEnsemble, getattr(ASMO, self.mdict["ensemble"]))
        )

        # add temperature if needed
        if self.mdict["temperature"] is not None:
            temperature = self.kg.create_node(
                f"temperature_{self.main_id}", ASMO.InputParameter
            )
            self.kg.add(
                (temperature, RDFS.label, Literal("temperature", datatype=XSD.string))
            )
            self.kg.add((activity, ASMO.hasInputParameter, temperature))
            self.kg.add(
                (
                    temperature,
                    ASMO.hasValue,
                    Literal(self.mdict["temperature"], datatype=XSD.float),
                )
            )
            self.kg.add(
                (temperature, ASMO.hasUnit, URIRef("http://qudt.org/vocab/unit/K"))
            )

        if self.mdict["pressure"] is not None:
            pressure = self.kg.create_node(
                f"pressure_{self.main_id}", ASMO.InputParameter
            )
            self.kg.add(
                (pressure, RDFS.label, Literal("pressure", datatype=XSD.string))
            )
            self.kg.add((activity, ASMO.hasInputParameter, pressure))
            self.kg.add(
                (
                    pressure,
                    ASMO.hasValue,
                    Literal(self.mdict["pressure"], datatype=XSD.float),
                )
            )
            self.kg.add(
                (pressure, ASMO.hasUnit, URIRef("http://qudt.org/vocab/unit/GigaPA"))
            )

        # potentials need to be mapped
        potential = URIRef(f"potential_{self.main_id}")
        if "meam" in self.mdict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.ModifiedEmbeddedAtomModel))
        elif "eam" in self.mdict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.EmbeddedAtomModel))
        elif "lj" in self.mdict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.LennardJonesPotential))
        elif "ace" in self.mdict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.MachineLearningPotential))
        else:
            self.kg.add((potential, RDF.type, ASMO.InteratomicPotential))

        if "uri" in self.mdict["potential"].keys():
            self.kg.add(
                (
                    potential,
                    CMSO.hasReference,
                    Literal(self.mdict["potential"]["uri"], datatype=XSD.string),
                )
            )
        if "label" in self.mdict["potential"].keys():
            self.kg.add(
                (potential, RDFS.label, Literal(self.mdict["potential"]["label"]))
            )

        self.kg.add((method, ASMO.hasInteratomicPotential, potential))
