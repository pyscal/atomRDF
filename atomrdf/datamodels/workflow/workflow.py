from typing import List, Optional, Union
import os
import numpy as np
import yaml
import uuid
import json
from pydantic import Field, field_validator
from atomrdf.datamodels.basemodels import (
    TemplateMixin,
    DataProperty,
    RDFMixin,
    BaseModel,
)
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.datamodels.activity import Activity
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
    MDO,
)
from atomrdf.datamodels.workflow.dof import *
from atomrdf.datamodels.workflow.ensemble import *
from atomrdf.datamodels.workflow.potential import *
from atomrdf.datamodels.workflow.software import *
from atomrdf.datamodels.workflow.method import *
from atomrdf.datamodels.workflow.algorithm import *
from atomrdf.datamodels.workflow.xcfunctional import *
from atomrdf.utils import get_simulation
from atomrdf.datamodels.workflow.property import *


class Simulation(Activity):
    label: Optional[str] = Field(default=None, description="Label of the simulation")
    # classes
    method: Optional[
        Union[
            MolecularDynamics,
            MolecularStatics,
            DensityFunctionalTheory,
        ]
    ] = Field(default=None, description="Computational method used in the simulation")

    algorithm: Optional[
        Union[
            EquationOfStateFit,
            QuasiHarmonicApproximation,
            ThermodynamicIntegration,
            ANNNIModel,
        ]
    ] = Field(default=None, description="Algorithm used in the simulation")

    # named individuals
    degrees_of_freedom: Optional[
        List[Union[AtomicPositionRelaxation, CellVolumeRelaxation, CellShapeRelaxation]]
    ] = Field(default=[], description="Degrees of freedom associated with the method")

    # named individuals
    thermodynamic_ensemble: Optional[
        Union[
            CanonicalEnsemble,
            MicrocanonicalEnsemble,
            IsothermalIsobaricEnsemble,
            IsoenthalpicIsobaricEnsemble,
            GrandCanonicalEnsemble,
        ]
    ] = Field(default=None, description="Thermodynamic ensemble used in the method")

    # class
    interatomic_potential: Optional[
        Union[
            InteratomicPotential,
            ModifiedEmbeddedAtomModel,
            EmbeddedAtomModel,
            LennardJonesPotential,
            MachineLearningPotential,
        ]
    ] = Field(default=None, description="Interatomic potential used in the method")

    # class
    xc_functional: Optional[Union[XCFunctional, GGA, LDA]] = Field(
        default=None, description="XC functional used in the method"
    )

    # class
    workflow_manager: Optional[SoftwareAgent] = Field(
        default=None, description="Workflow manager used in the simulation"
    )
    software: Optional[List[SoftwareAgent]] = Field(
        default=None, description="Softwares used in the simulation"
    )

    # list of classes
    input_parameter: Optional[List[InputParameter]] = Field(
        default=[], description="Input properties used in the simulation"
    )
    output_parameter: Optional[List[OutputParameter]] = Field(
        default=[], description="Output properties generated in the simulation"
    )
    calculated_property: Optional[List[CalculatedProperty]] = Field(
        default=[], description="Calculated properties from the simulation"
    )

    path: Optional[str] = Field(
        default=None, description="Path to the simulation directory"
    )

    @field_validator("method", mode="before")
    @classmethod
    def _validate_method(cls, v):
        if isinstance(v, str):
            if v in method_map:
                return method_map[v]()
            else:
                raise ValueError(f"Unknown method type: {v}")
        return v

    @field_validator("algorithm", mode="before")
    @classmethod
    def _validate_method(cls, v):
        if isinstance(v, str):
            if v in algorithm_map:
                return algorithm_map[v]()
            else:
                raise ValueError(f"Unknown algorithm type: {v}")
        return v

    @field_validator("degrees_of_freedom", mode="before")
    @classmethod
    def _validate_dof(cls, v):
        if isinstance(v, list):
            validated_dof = []
            for item in v:
                if isinstance(item, str):
                    if item in dof_map:
                        validated_dof.append(dof_map[item]())
                    else:
                        raise ValueError(f"Unknown degree of freedom: {item}")
                else:
                    raise ValueError(
                        f"Invalid type for degree of freedom: {type(item)}"
                    )
            return validated_dof
        return v

    @field_validator("thermodynamic_ensemble", mode="before")
    @classmethod
    def _validate_ensemble(cls, v):
        if isinstance(v, str):
            if v in ensemble_map:
                return ensemble_map[v]()
            else:
                raise ValueError(f"Unknown ensemble type: {v}")
        return v

    @field_validator("interatomic_potential", mode="before")
    @classmethod
    def _validate_potential(cls, v):
        if isinstance(v, dict):
            if v.get("potential_type") in potential_map:
                return potential_map[v.get("potential_type")](**v)
            else:
                return potential_map["InteratomicPotential"](**v)
        return v

    @field_validator("xc_functional", mode="before")
    @classmethod
    def _validate_xc_functional(cls, v):
        if isinstance(v, str):
            if v in xc_map:
                return xc_map[v]()
            else:
                return XCFunctional()
        return v

    @field_validator("software", mode="before")
    @classmethod
    def _validate_software(cls, v):
        if isinstance(v, list):
            return [SoftwareAgent(**item) for item in v]
        return [SoftwareAgent(**v)]

    @field_validator("workflow_manager", mode="before")
    @classmethod
    def _validate_workflow_manager(cls, v):
        if isinstance(v, dict):
            return SoftwareAgent(**v)
        return v

    def _to_graph_md_details(self, graph, simulation):
        # add ensemble
        if self.thermodynamic_ensemble:
            ensemble = self.thermodynamic_ensemble.to_graph(
                graph,
            )
            graph.add((simulation, ASMO.hasStatisticalEnsemble, ensemble))

        # add potential
        if self.interatomic_potential:
            potential = self.interatomic_potential.to_graph(graph, simulation)
            graph.add((simulation, ASMO.hasInteratomicPotential, potential))

    @classmethod
    def _from_graph_md_details(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        ensemble = graph.value(sim, ASMO.hasStatisticalEnsemble)
        if ensemble:
            cls.thermodynamic_ensemble = ensemble.toPython().split("/")[-1]

        potential = graph.value(sim, ASMO.hasInteratomicPotential)
        if potential:
            pot_type = potential.toPython().split("/")[-1]
            if pot_type == "ModifiedEmbeddedAtomModel":
                cls.interatomic_potential = ModifiedEmbeddedAtomModel.from_graph(
                    graph, potential
                )
            elif pot_type == "EmbeddedAtomModel":
                cls.interatomic_potential = EmbeddedAtomModel.from_graph(
                    graph, potential
                )
            elif pot_type == "LennardJonesPotential":
                cls.interatomic_potential = LennardJonesPotential.from_graph(
                    graph, potential
                )
            elif pot_type == "MachineLearningPotential":
                cls.interatomic_potential = MachineLearningPotential.from_graph(
                    graph, potential
                )
            else:
                cls.interatomic_potential = InteratomicPotential.from_graph(
                    graph, potential
                )
        return cls

    def _to_graph_dft_details(self, graph, simulation):
        # add XC functional
        if self.xc_functional:
            xc_functional = self.xc_functional.to_graph()
            graph.add((simulation, MDO.hasXCFunctional, xc_functional))

    @classmethod
    def _from_graph_dft_details(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        xc_functional = graph.value(sim, MDO.hasXCFunctional)
        if xc_functional:
            cls.xc_functional = xc_functional.toPython().split("/")[-1]
        return cls

    def _to_graph_dof(self, graph, simulation):
        if self.degrees_of_freedom:
            for dof in self.degrees_of_freedom:
                graph.add((simulation, ASMO.hasRelaxationDOF, dof.to_graph()))

    @classmethod
    def _from_graph_dof(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        dofs = [x[2] for x in graph.triples((sim, ASMO.hasRelaxationDOF, None))]
        doflist = []
        for dof in dofs:
            doflist.append(dof.toPython().split("/")[-1])

        if doflist:
            cls.degrees_of_freedom = doflist
        return cls

    def _to_graph_software(self, graph, simulation):
        workflow_manager = None
        if self.workflow_manager:
            workflow_manager = self.workflow_manager.to_graph(
                graph,
            )
            graph.add((simulation, PROV.wasAssociatedWith, workflow_manager))

        if self.software:
            for software in self.software:
                softobj = software.to_graph(
                    graph,
                )
                graph.add((simulation, PROV.wasAssociatedWith, softobj))
                if workflow_manager:
                    graph.add((workflow_manager, PROV.actedOnBehalfOf, softobj))

    @classmethod
    def _from_graph_software(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        workflow_manager = graph.value(sim, PROV.wasAssociatedWith)
        if workflow_manager:
            cls.workflow_manager = SoftwareAgent.from_graph(graph, workflow_manager)

        software = [x[2] for x in graph.triples((sim, PROV.wasAssociatedWith, None))]
        if software:
            cls.software = [SoftwareAgent.from_graph(graph, s) for s in software]
        return cls

    def to_graph_input_parameters(self, graph, simulation):
        if self.input_parameter:
            for param in self.input_parameter:
                param_uri = param.to_graph(graph)
                graph.add(
                    (simulation, ASMO.hasInputParameter, param_uri), validate=False
                )

    @classmethod
    def from_graph_input_parameters(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        input_params = [
            x[2] for x in graph.triples((sim, ASMO.hasInputParameter, None))
        ]
        if input_params:
            cls.input_parameter = [
                InputParameter.from_graph(graph, p) for p in input_params
            ]
        return cls

    def to_graph_output_parameters(self, graph, simulation):
        if self.output_parameter:
            for param in self.output_parameter:
                param_uri = param.to_graph(graph)
                graph.add(
                    (simulation, ASMO.hasOutputParameter, param_uri), validate=False
                )
                if param.associate_to_sample:
                    if self.output_sample:
                        output_samples = (
                            self.output_sample
                            if isinstance(self.output_sample, list)
                            else [self.output_sample]
                        )
                        for out_sample in output_samples:
                            graph.add(
                                (
                                    URIRef(out_sample),
                                    ASMO.hasCalculatedProperty,
                                    param_uri,
                                ),
                                validate=False,
                            )

    @classmethod
    def from_graph_output_parameters(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        output_params = [
            x[2] for x in graph.triples((sim, ASMO.hasOutputParameter, None))
        ]
        if output_params:
            cls.output_parameter = [
                OutputParameter.from_graph(graph, p) for p in output_params
            ]
        return cls

    @classmethod
    def from_graph_calculated_properties(cls, graph, sim_id):
        sim = get_simulation(graph, sim_id)
        calc_props = [x[2] for x in graph.triples((sim, ASMO.wasCalculatedBy, None))]
        if calc_props:
            cls.calculated_property = [
                CalculatedProperty.from_graph(graph, p) for p in calc_props
            ]
        return cls

    def to_graph_calculated_properties(self, graph, simulation):
        if self.calculated_property:
            for param in self.calculated_property:
                param_uri = param.to_graph(graph)
                graph.add((param_uri, ASMO.wasCalculatedBy, simulation), validate=False)
                if param.associate_to_sample:
                    if self.output_sample:
                        output_samples = (
                            self.output_sample
                            if isinstance(self.output_sample, list)
                            else [self.output_sample]
                        )
                        for out_sample in output_samples:
                            graph.add(
                                (
                                    URIRef(out_sample),
                                    ASMO.hasCalculatedProperty,
                                    param_uri,
                                ),
                                validate=False,
                            )

    def to_graph(self, graph):
        # if needed, serialise structures
        if self.input_sample:
            if isinstance(self.input_sample, list):
                self.input_sample = [
                    s.to_graph(graph) if isinstance(s, AtomicScaleSample) else s
                    for s in self.input_sample
                ]
            elif isinstance(self.input_sample, AtomicScaleSample):
                self.input_sample = self.input_sample.to_graph(graph)
        if self.output_sample:
            if isinstance(self.output_sample, list):
                self.output_sample = [
                    s.to_graph(graph) if isinstance(s, AtomicScaleSample) else s
                    for s in self.output_sample
                ]
            elif isinstance(self.output_sample, AtomicScaleSample):
                self.output_sample = self.output_sample.to_graph(graph)

        # create main simulation id
        main_id = uuid.uuid4()
        main_id = f"simulation:{main_id}"

        # add method
        method = self.method.to_graph(graph, main_id)

        # create simulation node based on method
        if self.method.basename in ["MolecularStatics", "MolecularDynamics"]:
            simulation = graph.create_node(main_id, ASMO.EnergyCalculation)
            graph.add((simulation, ASMO.hasComputationalMethod, method))
            self._to_graph_dof(graph, simulation)
            self._to_graph_md_details(graph, simulation)

        elif self.method.basename == "DensityFunctionalTheory":
            simulation = graph.create_node(main_id, ASMO.EnergyCalculation)
            graph.add((simulation, ASMO.hasComputationalMethod, method))
            self._to_graph_dof(graph, simulation)
            self._to_graph_dft_details(graph, simulation)

        else:
            simulation = graph.create_node(main_id, ASMO.Simulation)

        if self.algorithm.basename == "EquationOfStateFit":
            algorithm = graph.create_node(main_id, ASMO.EquationOfStateFit)

        elif self.algorithm.basename == "QuasiHarmonicApproximation":
            algorithm = graph.create_node(main_id, ASMO.QuasiHarmonicApproximation)

        elif self.method.basename == "ThermodynamicIntegration":
            algorithm = graph.create_node(main_id, ASMO.ThermodynamicIntegration)

        elif self.method.basename == "ANNNIModel":
            algorithm = graph.create_node(main_id, ASMO.ANNNImodel)

        graph.add((simulation, ASMO.usesSimulationAlgorithm, algorithm))

        # now add software
        self._to_graph_software(graph, simulation)

        # add structure layers
        if self.output_sample:
            output_samples = (
                self.output_sample
                if isinstance(self.output_sample, list)
                else [self.output_sample]
            )
            for out_sample in output_samples:
                graph.add((URIRef(out_sample), PROV.wasGeneratedBy, simulation))
                if self.input_sample:
                    input_samples = (
                        self.input_sample
                        if isinstance(self.input_sample, list)
                        else [self.input_sample]
                    )
                    for in_sample in input_samples:
                        graph.add(
                            (
                                URIRef(out_sample),
                                PROV.wasDerivedFrom,
                                URIRef(in_sample),
                            )
                        )

        if self.path:
            graph.add(
                (simulation, CMSO.hasPath, Literal(self.path, datatype=XSD.string))
            )

        if self.input_parameter:
            self.to_graph_input_parameters(graph, simulation)
        if self.output_parameter:
            self.to_graph_output_parameters(graph, simulation)
        if self.calculated_property:
            self.to_graph_calculated_properties(graph, simulation)

        return simulation

    @classmethod
    def from_graph(cls, graph, sim_id):
        cls = cls._from_graph_md_details(graph, sim_id)
        cls = cls._from_graph_dft_details(graph, sim_id)
        cls = cls._from_graph_dof(graph, sim_id)
        cls = cls._from_graph_software(graph, sim_id)
        cls = cls.from_graph_input_parameters(graph, sim_id)
        cls = cls.from_graph_output_parameters(graph, sim_id)
        cls = cls.from_graph_calculated_properties(graph, sim_id)

        sim = get_simulation(graph, sim_id)
        label = graph.get_label(sim)
        if label:
            cls.label = label

        # Get all input samples (can be multiple)
        input_samples = [
            str(x[2])
            for x in graph.triples((None, PROV.wasDerivedFrom, None))
            if any(graph.triples((x[0], PROV.wasGeneratedBy, sim)))
        ]
        if input_samples:
            cls.input_sample = (
                input_samples if len(input_samples) > 1 else input_samples[0]
            )

        # Get all output samples (can be multiple)
        output_samples = [
            str(x[0]) for x in graph.triples((None, PROV.wasGeneratedBy, sim))
        ]
        if output_samples:
            cls.output_sample = (
                output_samples if len(output_samples) > 1 else output_samples[0]
            )

        path = graph.value(sim, CMSO.hasPath)
        if path:
            cls.path = str(path)

        return cls
