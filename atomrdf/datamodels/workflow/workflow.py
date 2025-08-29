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
)
from atomrdf.datamodels.workflow.dof import *
from atomrdf.datamodels.workflow.ensemble import *
from atomrdf.datamodels.workflow.potential import *
from atomrdf.datamodels.workflow.software import *
from atomrdf.datamodels.workflow.method import *


class Simulation(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the method")
    label: Optional[str] = Field(default=None, description="Label of the method")
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

    degrees_of_freedom: Optional[
        List[Union[AtomicPositionRelaxation, CellVolumeRelaxation, CellShapeRelaxation]]
    ] = Field(default=[], description="Degrees of freedom associated with the method")

    thermodynamic_ensemble: Optional[
        List[
            Union[
                CanonicalEnsemble,
                MicrocanonicalEnsemble,
                IsothermalIsobaricEnsemble,
                IsoenthalpicIsobaricEnsemble,
                GrandCanonicalEnsemble,
            ]
        ]
    ] = Field(default=None, description="Thermodynamic ensemble used in the method")

    interatomic_potential: Optional[
        Union[
            InteratomicPotential,
            ModifiedEmbeddedAtomModel,
            EmbeddedAtomModel,
            LennardJonesPotential,
            MachineLearningPotential,
        ]
    ] = Field(default=None, description="Interatomic potential used in the method")

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
        if isinstance(v, list):
            ensemble_map = {
                "CanonicalEnsemble": CanonicalEnsemble,
                "MicrocanonicalEnsemble": MicrocanonicalEnsemble,
                "IsothermalIsobaricEnsemble": IsothermalIsobaricEnsemble,
                "IsoenthalpicIsobaricEnsemble": IsoenthalpicIsobaricEnsemble,
                "GrandCanonicalEnsemble": GrandCanonicalEnsemble,
            }
            validated_ensembles = []
            for item in v:
                if isinstance(item, str):
                    if item in ensemble_map:
                        validated_ensembles.append(ensemble_map[item]())
                    else:
                        raise ValueError(f"Unknown ensemble type: {item}")
                else:
                    raise ValueError(f"Invalid type for ensemble: {type(item)}")
            return validated_ensembles
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

    def _add_md_details(self, graph, simulation):
        if self.thermodynamic_ensemble:
            for ensemble in self.thermodynamic_ensemble:
                graph.add(
                    (
                        simulation,
                        ASMO.hasStatisticalEnsemble,
                        getattr(ASMO, ensemble.basename),
                    )
                )

    def _add_dof(self, graph, simulation):
        if self.degrees_of_freedom:
            for dof in self.degrees_of_freedom:
                graph.add(
                    (simulation, ASMO.hasRelaxationDOF, getattr(ASMO, dof.basename))
                )

    def to_graph(self, graph, main_id):
        simulation = graph.create_node(main_id, ASMO.EnergyCalculation)
        graph.add((simulation, ASMO.hasComputationalMethod, ASMO.MolecularStatics))

        # add DOF
        self._add_dof(graph, simulation)
        self._add_md_details(graph, simulation)

        if self.interatomic_potential:
            potential = self.interatomic_potential.to_graph(graph, main_id)
            graph.add(
                (
                    simulation,
                    ASMO.hasInteratomicPotential,
                    potential,
                )
            )
