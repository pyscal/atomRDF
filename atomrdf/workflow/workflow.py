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

from rdflib import Namespace, XSD, RDF, RDFS, BNode, URIRef

import warnings
import numpy as np
import os
import copy
import ast
import uuid
import importlib
import uuid

from atomrdf.structure import System

# Move imports to another file
from atomrdf.namespace import PROV, CMSO, PODO, ASMO, MDO, Literal, UNSAFECMSO, UNSAFEASMO

class Workflow:
    def __init__(self, kg):
        """
        Initialize the workflow environment

        Parameters
        ----------
        kg: pyscal-rdf KnowledgeGraph
        environment: string
            the workflow environment. This is used to import the necessary functions.

        """
        self.kg = kg

    def inform_graph(self, pr, workflow_environment=None, workflow_module=None):
        if workflow_environment is not None:
            workflow_module = importlib.import_module(f"atomrdf.workflow.{workflow_environment}.{workflow_environment}")
        
        if workflow_module is not None:
            workflow_module.inform_graph(pr, self.kg)
        

    def to_graph(self, job=None, workflow_environment=None, workflow_module=None, job_dicts=None, add_intermediate_jobs=False):

        if workflow_environment is not None:
            workflow_module = importlib.import_module(f"atomrdf.workflow.{workflow_environment}.{workflow_environment}")
            job_dicts = np.atleast_1d(workflow_module.process_job(job))
        elif workflow_module is not None:
            job_dicts = np.atleast_1d(workflow_module.process_job(job))
        if job_dicts is None:
            raise ValueError("Job dict could not be calculated!")
        job_dicts = np.atleast_1d(job_dicts)
        
        #print(job_dict)

        #now we call the functions in order
        for job_dict in job_dicts:
            if not 'intermediate' in job_dict.keys():
                job_dict['intermediate'] = False
            if (not add_intermediate_jobs) and job_dict['intermediate']:
                continue
            job_dict = self._add_structure(job_dict)
            parent_sample = job_dict['sample']['initial']
            sample = job_dict['sample']['final']
            structure = job_dict['structure']['final']
            #print(parent_sample)
            self._add_structural_relation(parent_sample, sample, structure)
            self._add_method(job_dict)

        #now add child connection
        if add_intermediate_jobs:
            for job_dict in job_dicts:
                if not job_dict['intermediate']:
                    parent = job_dict
                    break
            
            for job_dict in job_dicts:
                if job_dict['intermediate']:
                    parent_sample = job_dict['sample']['final']
                    sample = parent['sample']['final']
                    structure = parent['structure']['final']
                    self._add_structural_relation(parent_sample, sample, structure)

    def _add_structure(self, job_dict):
        #ensure these are not strings
        if isinstance(job_dict['sample']['initial'], str):
            job_dict['sample']['initial'] = URIRef(job_dict['sample']['initial'])
        if isinstance(job_dict['sample']['final'], str):
            job_dict['sample']['final'] = URIRef(job_dict['sample']['final'])

        if (job_dict['structure']['final'] is None) and (job_dict['sample']['final'] is None):
            raise ValueError("Either structure or sample should be specified")

        if job_dict['sample']['final'] is None:
            # its not added to graph yet
            structure = job_dict['structure']['final']
            structure.graph = self.kg
            structure.to_graph()
            job_dict['sample']['final'] = structure.sample


        if (job_dict['sample']['initial'] is None) or (job_dict['sample']['initial'] not in self.kg.sample_ids):
            # its not added to graph yet
            parent_structure = job_dict['structure']['initial']
            if parent_structure is not None:
                parent_structure.graph = self.kg
                parent_structure.to_graph()
                job_dict['sample']['initial'] = parent_structure.sample
        return job_dict

    def _add_structural_relation(
        self, parent_sample, sample, structure,
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

        self.kg.add((sample, RDF.type, PROV.Entity))
        
        if parent_sample is not None:
            self.kg.add((parent_sample, RDF.type, PROV.Entity))
            self.kg.add((sample, PROV.wasDerivedFrom, parent_sample))
            #update label
            label = self.kg.get_string_label(sample)
            parent_label = self.kg.get_string_label(parent_sample)
            if label is None:
                new_label = f"{parent_label}_derivative"
            else:
                new_label = f"{label}_from_{parent_label}"
            self.kg.change_label(sample, new_label)
            self._get_lattice_properties(parent_sample, sample, structure,)
            self._add_inherited_properties(parent_sample, sample,)
            self._add_cell_repetitions(parent_sample, sample,)

    def _get_lattice_properties(
        self, parent_sample, sample, structure,
    ):
        if sample is None:
            return

        material = list(
            [
                k[2]
                for k in self.kg.triples((parent_sample, CMSO.hasMaterial, None))
            ]
        )[0]
        crystal_structure = self.kg.value(material, CMSO.hasStructure)

        altname = self.kg.value(crystal_structure, CMSO.hasAltName)

        #space_group = self.kg.value(crystal_structure, CMSO.hasSpaceGroup)
        space_group_symbol = self.kg.value(crystal_structure, CMSO.hasSpaceGroupSymbol)
        space_group_number = self.kg.value(crystal_structure, CMSO.hasSpaceGroupNumber)

        unit_cell = self.kg.value(crystal_structure, CMSO.hasUnitCell)
        blattice = self.kg.value(
            unit_cell,
            CMSO.hasBravaisLattice,
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

        structure._add_crystal_structure(targets=targets)

    def _add_cell_repetitions(self, parent_sample, sample):
        #add lattice repetitions when doing simulations

        if sample is None:
            return
        
        old_simcell = self.kg.value(parent_sample, CMSO.hasSimulationCell)
        x = self.kg.value(old_simcell, CMSO.hasRepetition_x)        
        y = self.kg.value(old_simcell, CMSO.hasRepetition_y)
        z = self.kg.value(old_simcell, CMSO.hasRepetition_z)
        x = 1 if x is None else x
        y = 1 if y is None else y
        z = 1 if z is None else z

        new_simcell = self.kg.value(sample, CMSO.hasSimulationCell)
        self.kg.add((new_simcell, CMSO.hasRepetition_x, x))
        self.kg.add((new_simcell, CMSO.hasRepetition_y, y))
        self.kg.add((new_simcell, CMSO.hasRepetition_z, z))
    
    def _add_inherited_properties(
        self, parent_sample, sample,
    ):
        # Here we need to add inherited info: CalculatedProperties will be lost
        # Defects will be inherited

        if sample is None:
            return

        self.kg.copy_defects(sample, parent_sample)

    def _add_method(
            self, job_dict, 
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

        # add activity
        # ----------------------------------------------------------
        main_id = uuid.uuid4()
        main_id = f'activity:{main_id}'
        job_dict['id'] = main_id
        activity = URIRef(main_id)
        self.kg.add((activity, RDF.type, PROV.Activity))
        
        # add method
        # ----------------------------------------------------------
        method = URIRef(f"{main_id}_method")
        if job_dict["method"] == "MolecularStatics":
            self.kg.add((activity, RDF.type, ASMO.EnergyCalculation))
            self.kg.add((method, RDF.type, ASMO.MolecularStatics))
            self._add_dof(job_dict, activity)
            self._add_md(job_dict, activity)

        elif job_dict["method"] == "MolecularDynamics":
            self.kg.add((activity, RDF.type, ASMO.EnergyCalculation))
            self.kg.add((method, RDF.type, ASMO.MolecularDynamics))
            self._add_dof(job_dict, activity)
            self._add_md(job_dict, activity)

        elif job_dict["method"] == "DensityFunctionalTheory":
            self.kg.add((activity, RDF.type, ASMO.EnergyCalculation))
            self.kg.add((method, RDF.type, ASMO.DensityFunctionalTheory))
            self._add_dof(job_dict, activity)
            self._add_dft(job_dict, method, activity)

        elif job_dict["method"] == "EquationOfState":            
            #special type of EOS should be initialised!
            self.kg.add((activity, RDF.type, ASMO.EquationOfStateFit))

        elif job_dict["method"] == "QuasiHarmonicModel":            
            self.kg.add((activity, RDF.type, ASMO.QuasiHarmonicModel))

        elif job_dict["method"] == "ThermodynamicIntegration":     
            self._add_dof(job_dict, activity)
            self._add_md(job_dict, activity)       
            self.kg.add((activity, RDF.type, ASMO.ThermodynamicIntegration))

        # add that structure was generated
        self.kg.add((activity, ASMO.hasComputationalMethod, method))
        self.kg.add((job_dict['sample']['final'], PROV.wasGeneratedBy, activity))
        if 'path' in job_dict.keys():
            self.kg.add((activity, CMSO.hasPath, Literal(job_dict['path'], datatype=XSD.string)))
        self._add_inputs(job_dict, activity)
        self._add_outputs(job_dict, activity)
        self._add_software(job_dict, method)

    def _add_dof(self, job_dict, activity):
        for dof in job_dict["dof"]:
            self.kg.add((activity, ASMO.hasRelaxationDOF, getattr(ASMO, dof)))
    
    def _select_base_property(self, out, main_id, default_class):
        if "base" in out.keys():
            base = out["base"]
        else:
            base = out["label"]

        if base == 'TotalEnergy':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.TotalEnergy,
                label=out["label"],
            )
        elif base == 'PotentialEnergy':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.PotentialEnergy,
                label=out["label"],
            )
        elif base == 'KineticEnergy':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.KineticEnergy,
                label=out["label"],
            )
        elif base == 'Volume':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.Volume,
                label=out["label"],
            )
        elif base == 'Pressure':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.Pressure,
                label=out["label"],
            )
        elif base == 'Temperature':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.Temperature,
                label=out["label"],
            )
        elif base == 'BulkModulus':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.BulkModulus,
                label=out["label"],
            )
        elif base == 'FreeEnergy':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.FreeEnergy,
                label=out["label"],
            )
        elif base == 'EnergyCutoff':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.EnergyCutoff,
                label=out["label"],
            )
        elif base == 'ExplicitKPointMesh':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.ExplicitKPointMesh,
                label=out["label"],
            )
        elif base == 'GammaCenteredKPointMesh':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.GammaCenteredKPointMesh,
                label=out["label"],
            )
        elif base == 'MonkhorstPackKPointMesh':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.MonkhorstPackKPointMesh,
                label=out["label"],
            )
        elif base == 'KPointMesh':
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', UNSAFEASMO.KPointMesh,
                label=out["label"],
            )
        else:
            prop = self.kg.create_node(
                f'{main_id}_{out["label"]}', default_class,
                label=out["label"],
            )
        
        self.kg.add((prop, UNSAFEASMO.hasValue, Literal(out["value"])))
        if "unit" in out.keys():
            unit = out["unit"]
            self.kg.add(
                (
                    prop,
                    UNSAFEASMO.hasUnit,
                    URIRef(f"http://qudt.org/vocab/unit/{unit}"),
                )
            )
        return prop


    def _add_inputs(self, job_dict, activity):
        main_id = job_dict['id']
        if "inputs" in job_dict.keys():
            for inp in job_dict["inputs"]:
                prop = self._select_base_property(inp, main_id, ASMO.InputParameter)
                self.kg.add((activity, UNSAFEASMO.hasInputParameter, prop))

    def _add_outputs(self, job_dict, activity):
        main_id = job_dict['id']
        if "outputs" in job_dict.keys():
            for out in job_dict["outputs"]:
                #here we add the classes by property
                #call func here
                prop = self._select_base_property(out, main_id, CMSO.CalculatedProperty)
                self.kg.add((prop, UNSAFEASMO.wasCalculatedBy, activity))
                
                if out["associate_to_sample"]:
                    self.kg.add((job_dict['sample']['final'], UNSAFECMSO.hasCalculatedProperty, prop))

    def _add_software(self, job_dict, method):
        # finally add software
        wfagent = None
        if "workflow_manager" in job_dict.keys():
            wfagent = self.kg.create_node(
                job_dict["workflow_manager"]["uri"], PROV.SoftwareAgent
            )
            self.kg.add(
                (wfagent, RDFS.label, Literal(job_dict["workflow_manager"]["label"]))
            )
            self.kg.add((method, PROV.wasAssociatedWith, wfagent))

        for software in job_dict["software"]:
            agent = self.kg.create_node(software["uri"], PROV.SoftwareAgent)
            self.kg.add((agent, RDFS.label, Literal(software["label"])))
            if wfagent is not None:
                self.kg.add((wfagent, PROV.actedOnBehalfOf, agent))
            else:
                self.kg.add((method, PROV.wasAssociatedWith, agent))


    def _add_dft(self, job_dict, method, activity):
        main_id = job_dict['id']
        if job_dict["xc_functional"] is not None:
            if job_dict["xc_functional"] in ['PBE', 'GGA']:
                self.kg.add((method, MDO.hasXCFunctional, MDO.GGA))
            elif job_dict["xc_functional"] in ['LDA']:
                self.kg.add((method, MDO.hasXCFunctional, MDO.LDA))

    
    def _add_md(self, job_dict, activity):
        main_id = job_dict['id']
        if job_dict["ensemble"] is not None:
            self.kg.add(
                (activity, ASMO.hasStatisticalEnsemble, getattr(ASMO, job_dict["ensemble"]))
            )

        # potentials need to be mapped
        potential = URIRef(f"{main_id}_potential")
        if "meam" in job_dict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.ModifiedEmbeddedAtomModel))
        elif "eam" in job_dict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.EmbeddedAtomModel))
        elif "lj" in job_dict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.LennardJonesPotential))
        elif "ace" in job_dict["potential"]["type"]:
            self.kg.add((potential, RDF.type, ASMO.MachineLearningPotential))
        else:
            self.kg.add((potential, RDF.type, ASMO.InteratomicPotential))

        if "uri" in job_dict["potential"].keys():
            self.kg.add(
                (
                    potential,
                    CMSO.hasReference,
                    Literal(job_dict["potential"]["uri"], datatype=XSD.string),
                )
            )
        if "label" in job_dict["potential"].keys():
            self.kg.add(
                (potential, RDFS.label, Literal(job_dict["potential"]["label"]))
            )

        self.kg.add((activity, ASMO.hasInteratomicPotential, potential))