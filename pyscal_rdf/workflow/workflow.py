"""
Workflows aspects for non-automated annotation of structures.

"""

from pyscal_rdf.structure import System
from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

import pyscal_rdf.workflow.pyiron as pi

import warnings
import numpy as np
import os
import copy
import ast
import uuid

PROV = Namespace("http://www.w3.org/ns/prov#")
CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")
ASO = Namespace("http://purls.helmholtz-metadaten.de/aso/")

class Workflow:
    def __init__(self, kg, 
        environment='pyiron'):
        self.kg = kg
        if environment == 'pyiron':
            self.wenv = pi
        else:
            raise ValueError('unknow workflow environment')

    def _prepare_job(self, workflow_object):
        self.wenv._check_if_job_is_valid(workflow_object)
        parent_structure, parent_sample, structure, sample = self.wenv._add_structures(self.kg, workflow_object)
        method_dict = self.wenv._identify_method(workflow_object)

        if (structure is None) and (sample is None):
            raise ValueError('Either structure or sample should be specified')

        if sample is None:
            #its not added to graph yet
            structure.graph = self.kg
            structure.to_graph()
            sample = structure.sample
        
        if parent_sample is None:
            #its not added to graph yet
            parent_structure.graph = self.kg
            parent_structure.to_graph()
            parent_sample = parent_structure.sample

        self.sample = sample
        self.mdict = method_dict
        self.parent_sample = parent_sample

    def _add_inherited_properties(self, ):
        #Here we need to add inherited info: CalculatedProperties will be lost
        #Defects will be inherited
        if self.parent_sample is None:
            return

        parent_material = list([k[2] for k in self.kg.graph.triples((self.parent_sample, CMSO.hasMaterial, None))])[0]
        parent_defects = list([x[2] for x in self.kg.graph.triples((parent_material, CMSO.hasDefect, None))])
        #now for each defect we copy add this to the final sample
        material = list([k[2] for k in self.kg.graph.triples((self.sample, CMSO.hasMaterial, None))])[0]

        for defect in parent_defects:
            new_defect = URIRef(defect.toPython())
            self.kg.graph.add((material, CMSO.hasDefect, new_defect))
            #now fetch all defect based info
            for triple in self.kg.graph.triples((defect, None, None)):
                self.kg.graph.add((new_defect, triple[1], triple[2]))

        #now add the special props for vacancy
        parent_simcell = self.kg.graph.value(self.sample, CMSO.hasSimulationCell)
        simcell = self.kg.graph.value(self.parent_sample, CMSO.hasSimulationCell) 
        
        for triple in self.kg.graph.triples((parent_simcell, PODO.hasVacancyConcentration, None)):
            self.kg.graph.add((simcell, triple[1], triple[2]))
        for triple in self.kg.graph.triples((parent_simcell, PODO.hasNumberOfVacancies, None)):
            self.kg.graph.add((simcell, triple[1], triple[2]))

    def _get_lattice_properties(self, ):
        if self.parent_sample is None:
            return

        parent_material = list([k[2] for k in self.kg.graph.triples((self.parent_sample, CMSO.hasMaterial, None))])[0]
        parent_crystal_structure = self.kg.graph.value(parent_material, CMSO.hasStructure)
        parent_altname = self.kg.graph.value(parent_crystal_structure, CMSO.hasAltName)

        #add this to new structure
        material = list([k[2] for k in self.kg.graph.triples((self.sample, CMSO.hasMaterial, None))])[0]
        crystal_structure = self.kg.graph.value(material, CMSO.hasStructure)
        self.kg.add((crystal_structure, CMSO.hasAltName, parent_altname))

        #space group
        parent_space_group = self.kg.graph.value(parent_crystal_structure, CMSO.hasSpaceGroup)
        space_group = self.kg.graph.value(crystal_structure, CMSO.hasSpaceGroup)
        for triple in self.kg.graph.triples((parent_space_group, None, None)):
            self.kg.graph.add((space_group, triple[1], triple[2]))

        #unit cell
        parent_unit_cell = self.kg.graph.value(parent_crystal_structure, CMSO.hasUnitCell)
        parent_bv = self.kg.graph.value(parent_unit_cell, CMSO.hasBravaisLattice)

        unit_cell = self.kg.graph.value(crystal_structure, CMSO.hasUnitCell)
        self.kg.graph.add((unit_cell, CMSO.hasBravaisLattice, parent_bv))

        #lattice parameter
        parent_lattice_parameter = self.kg.graph.value(parent_unit_cell, CMSO.hasLatticeParameter)
        lattice_parameter = self.kg.graph.value(unit_cell, CMSO.hasLatticeParameter)
        for triple in self.kg.graph.triples((parent_lattice_parameter, None, None)):
            self.kg.graph.add((lattice_parameter, triple[1], triple[2]))

        #lattice angle
        parent_lattice_angle = self.kg.graph.value(parent_unit_cell, CMSO.hasAngle)
        lattice_angle = self.kg.graph.value(unit_cell, CMSO.hasAngle)
        for triple in self.kg.graph.triples((parent_lattice_angle, None, None)):
            self.kg.graph.add((lattice_angle, triple[1], triple[2]))


    def add_structural_relation(self, ):
        self.kg.add((self.sample, RDF.type, PROV.Entity))
        if self.parent_sample is not None:
            self.kg.add((self.parent_sample, RDF.type, PROV.Entity))
            self.kg.add((self.sample, PROV.wasDerivedFrom, self.parent_sample))
            self._add_inherited_properties()
            self._get_lattice_properties()


    def add_method(self, ):
        """
        mdict
        -----
        md:
           method: MolecularStatics
           temperature: 100
           pressure: 0
           dof:
             - AtomicPositions
             - CellVolume
           ensemble: NPT
           id: 2314
           potential:
             uri: https://doi.org/xxx
             type: eam
             label: string
           workflow_manager:
             uri: xxxx
             label: pyiron
           software:
           - uri: xxxx
             label: lammps
           - uri: xxxx
             label: pyscal

        """
        if self.mdict is None:
            return

        if 'md' in self.mdict.keys():
            method_type = 'md'
            mdict = self.mdict['md']
        elif 'dft' in self.mdict.keys():
            method_type = 'dft'
            mdict = self.mdict['dft']
        else:
            raise KeyError('method dict keys should be either md or dft')

        
        #add activity
        main_id = mdict['id']
        activity = URIRef(f'activity:{main_id}')
        self.kg.add((activity, RDF.type, PROV.Activity))

        #method, this is specific to dft/md
        if method_type == 'md':
            method = URIRef(f'method:{main_id}')
            if mdict['method'] == 'MolecularStatics':
                self.kg.add((method, RDF.type, ASO.MolecularStatics))
            elif mdict['method'] == 'MolecularDynamics':
                self.kg.add((method, RDF.type, ASO.MolecularDynamics))
        elif method_type == 'dft':
            method = URIRef(f'method:{main_id}')
            if mdict['method'] == 'DensityFunctionalTheory':
                self.kg.add((method, RDF.type, ASO.DensityFunctionalTheory))
        self.kg.add((activity, ASO.hasMethod, method))

        if len(mdict['dof']) == 0:
            self.kg.add((activity, RDF.type, ASO.RigidEnergyCalculation))
        else:
            self.kg.add((activity, RDF.type, ASO.StructureOptimization))

        for dof in mdict['dof']:
            self.kg.add((activity, ASO.hasRelaxationDOF, getattr(ASO, dof)))

        if method_type == 'md':
            self.kg.add((method, ASO.hasStatisticalEnsemble, getattr(ASO, mdict['ensemble'])))

            #add temperature if needed
            if mdict['temperature'] is not None:
                temperature = URIRef(f'temperature:{main_id}')
                self.kg.add((temperature, RDF.type, ASO.InputParameter))
                self.kg.add((temperature, RDFS.label, Literal('temperature', datatype=XSD.string)))
                self.kg.add((activity, ASO.hasInputParameter, temperature))
                self.kg.add((temperature, ASO.hasValue, Literal(mdict['temperature'], datatype=XSD.float)))
                self.kg.add((temperature, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))

            if mdict['pressure'] is not None:
                pressure = URIRef(f'pressure:{main_id}')
                self.kg.add((pressure, RDF.type, ASO.InputParameter))
                self.kg.add((pressure, RDFS.label, Literal('pressure', datatype=XSD.string)))
                self.kg.add((activity, ASO.hasInputParameter, pressure))
                self.kg.add((pressure, ASO.hasValue, Literal(mdict['pressure'], datatype=XSD.float)))
                self.kg.add((pressure, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/GigaPA')))

            #potentials need to be mapped
            potential = URIRef(f'potential:{main_id}')
            if 'meam' in mdict['potential']['type']:
                self.kg.add((potential, RDF.type, ASO.MEAM))
            elif 'eam' in mdict['potential']['type']:
                self.kg.add((potential, RDF.type, ASO.EAM))
            elif 'lj' in mdict['potential']['type']:
                self.kg.add((potential, RDF.type, ASO.LennardJones))
            elif 'ace' in mdict['potential']['type']:
                self.kg.add((potential, RDF.type, ASO.MLPotential))
            else:
                self.kg.add((potential, RDF.type, ASO.InteratomicPotential))

            if 'uri' in mdict['potential'].keys():
                self.kg.add((potential, ASO.hasReference, Literal(mdict['potential']['uri'])))
            if 'label' in mdict['potential'].keys():
                self.kg.add((potential, RDFS.label, Literal(mdict['potential']['label'])))

            self.kg.add((method, ASO.hasInteratomicPotential, potential))

        self.kg.add((self.sample, PROV.wasGeneratedBy, activity))

        #finally add software
        wfagent = None
        if 'workflow_manager' in mdict.keys():
            wfagent = URIRef(mdict["workflow_manager"]['uri'])
            self.kg.add((wfagent, RDF.type, PROV.SoftwareAgent))
            self.kg.add((wfagent, RDFS.label, Literal(mdict["workflow_manager"]['label'])))
            self.kg.add((method, PROV.wasAssociatedWith, wfagent))
        
        for software in mdict['software']:
            agent = URIRef(software['uri'])
            self.kg.add((agent, RDF.type, PROV.SoftwareAgent))
            self.kg.add((agent, RDFS.label, Literal(software['label'])))

            if wfagent is not None:
                self.kg.add((wfagent, PROV.actedOnBehalfOf, agent))
            else:
                self.kg.add((method, PROV.wasAssociatedWith, agent))

    def to_graph(self, workflow_object):
        self._prepare_job(workflow_object)
        self.add_structural_relation()
        self.add_method()