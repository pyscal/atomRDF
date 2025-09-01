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
    Activity,
)
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
from atomrdf.datamodels.workflow.xcfunctional import *


class Simulation(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the method")
    label: Optional[str] = Field(default=None, description="Label of the method")
    # classes
    method: Optional[
        Union[
            MolecularDynamics,
            MolecularStatics,
            DensityFunctionalTheory,
            EquationOfStateFit,
            QuasiHarmonicApproximation,
            ThermodynamicIntegration,
        ]
    ] = Field(default=None, description="Computational method used in the simulation")

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

    @field_validator("method", mode="before")
    @classmethod
    def _validate_method(cls, v):
        if isinstance(v, str):
            method_map = {
                "MolecularDynamics": MolecularDynamics,
                "MolecularStatics": MolecularStatics,
                "DensityFunctionalTheory": DensityFunctionalTheory,
                "EquationOfState": EquationOfStateFit,
                "QuasiHarmonicModel": QuasiHarmonicApproximation,
                "ThermodynamicIntegration": ThermodynamicIntegration,
            }
            if v in method_map:
                return method_map[v]()
            else:
                raise ValueError(f"Unknown method type: {v}")
        return v

    @field_validator("degrees_of_freedom", mode="before")
    @classmethod
    def _validate_dof(cls, v):
        if isinstance(v, list):
            dof_map = {
                "AtomicPositionRelaxation": AtomicPositionRelaxation,
                "CellVolumeRelaxation": CellVolumeRelaxation,
                "CellShapeRelaxation": CellShapeRelaxation,
            }
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
            ensemble_map = {
                "CanonicalEnsemble": CanonicalEnsemble,
                "MicrocanonicalEnsemble": MicrocanonicalEnsemble,
                "IsothermalIsobaricEnsemble": IsothermalIsobaricEnsemble,
                "IsoenthalpicIsobaricEnsemble": IsoenthalpicIsobaricEnsemble,
                "GrandCanonicalEnsemble": GrandCanonicalEnsemble,
            }
            if v in ensemble_map:
                return ensemble_map[v]()
            else:
                raise ValueError(f"Unknown ensemble type: {v}")
        return v

    @field_validator("interatomic_potential", mode="before")
    @classmethod
    def _validate_potential(cls, v):
        if isinstance(v, dict):
            potential_map = {
                "InteratomicPotential": InteratomicPotential,
                "ModifiedEmbeddedAtomModel": ModifiedEmbeddedAtomModel,
                "EmbeddedAtomModel": EmbeddedAtomModel,
                "LennardJonesPotential": LennardJonesPotential,
                "MachineLearningPotential": MachineLearningPotential,
            }
            if v.get("potential_type") in potential_map:
                return potential_map[v.get("potential_type")](**v)
            else:
                return potential_map["InteratomicPotential"](**v)
        return v

    @field_validator("xc_functional", mode="before")
    @classmethod
    def _validate_xc_functional(cls, v):
        if isinstance(v, str):
            xc_map = {
                "LDA": LDA,
                "GGA": GGA,
                "PBE": GGA,
            }
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

    def _add_md_details(self, graph, simulation):
        # add ensemble
        if self.thermodynamic_ensemble:
            ensemble = self.thermodynamic_ensemble.to_graph(graph, simulation)
            graph.add((simulation, ASMO.hasStatisticalEnsemble, ensemble))

        # add potential
        if self.interatomic_potential:
            potential = self.interatomic_potential.to_graph(graph, simulation)
            graph.add((simulation, ASMO.hasInteratomicPotential, potential))

    def _add_dft_details(self, graph, simulation):
        # add XC functional
        if self.xc_functional:
            xc_functional = self.xc_functional.to_graph()
            graph.add((simulation, MDO.hasXCFunctional, xc_functional))

    def _add_dof(self, graph, simulation):
        if self.degrees_of_freedom:
            for dof in self.degrees_of_freedom:
                graph.add((simulation, ASMO.hasRelaxationDOF, dof.to_graph()))

    def _add_software(self, graph, simulation):
        if self.software:
            software = self.software.to_graph()
            graph.add((simulation, ASMO.hasSoftware, software))
