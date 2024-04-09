"""
Workflows aspects for non-automated annotation of structures.

"""

from pyscal_rdf.structure import System
from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

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


def annotate_md(graph,
                structure,
                id = None,
                pressure = None,
                ensemble = None,
                temperature = None,
                potential_doi = None,
                potential_type = "",
                software_id = None,
                software = None,
                ):
    """
    Annotate a given structure with MD simulation details
    """
    main_id = str(uuid.uuid4())
    if id is not None:
        main_id = id
        
    sample = structure.sample
    graph.add((sample, RDF.type, PROV.Entity))

    activity = URIRef(f'activity:{main_id}')
    graph.add((activity, RDF.type, PROV.Activity))
    graph.add((activity, RDF.type, ASO.StructureOptimization))

    method = URIRef(f'method:{main_id}')
    graph.add((method, RDF.type, ASO.MolecularDynamics))
    graph.add((activity, ASO.hasMethod, method))
    graph.add((activity, ASO.hasRelaxationDOF, ASO.AtomicPosition))

    if pressure is None:
        pass
    elif np.isscalar(pressure):
        graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
    else: 
        #check if pressure is hydrostatic or not
        axial_all_alike = None not in pressure[:3] and np.allclose(
            pressure[:3], pressure[0]
        )
        shear_all_none = all(p is None for p in pressure[3:])
        shear_all_zero = None not in pressure[3:] and np.allclose(pressure[3:], 0)
        hydrostatic = axial_all_alike and (shear_all_none or shear_all_zero)
        if hydrostatic:
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
        else:
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellShape))

    if ensemble is not None:
        if ensemble == 'NVT':
            graph.add((method, ASO.hasStatisticalEnsemble, ASO.NVT))
            if temperature is None:
                raise ValueError('Temperature cannot be None in NVT')
            temperature_node = URIRef(f'temperature:{main_id}')
            graph.add((temperature_node, RDF.type, ASO.InputParameter))
            graph.add((temperature_node, RDFS.label, Literal('temperature', datatype=XSD.string)))
            graph.add((activity, ASO.hasInputParameter, temperature_node))
            graph.add((temperature_node, ASO.hasValue, Literal(temperature, datatype=XSD.float)))
            graph.add((temperature_node, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))

        elif ensemble == 'NPT':
            graph.add((method, ASO.hasStatisticalEnsemble, ASO.NPT))
            if temperature is None:
                raise ValueError('Temperature cannot be None in NPT')
            temperature_node = URIRef(f'temperature:{main_id}')
            graph.add((temperature_node, RDF.type, ASO.InputParameter))
            graph.add((temperature_node, RDFS.label, Literal('temperature', datatype=XSD.string)))
            graph.add((activity, ASO.hasInputParameter, temperature_node))
            graph.add((temperature_node, ASO.hasValue, Literal(temperature, datatype=XSD.float)))
            graph.add((temperature_node, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))

            pressure_node = URIRef(f'pressure:{main_id}')
            graph.add((pressure_node, RDF.type, ASO.InputParameter))
            graph.add((pressure_node, RDFS.label, Literal('pressure', datatype=XSD.string)))
            graph.add((activity, ASO.hasInputParameter, pressure_node))
            graph.add((pressure_node, ASO.hasValue, Literal(pressure, datatype=XSD.float)))
            graph.add((pressure_node, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/GigaPA')))

    graph.add((sample, PROV.wasGeneratedBy, activity))

    if potential_doi is None:
        warnings.warn('potential_doi is None, maybe consider providing it?')
    else:
        potential = URIRef(f'potential:{main_id}')

        if 'meam' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.MEAM))
        elif 'eam' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.EAM))
        elif 'lj' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.LennardJones))
        elif 'ace' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.MLPotential))
        elif 'snap' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.MLPotential))
        elif 'tersoff' in potential_type.lower():
            graph.add((potential, RDF.type, ASO.Tersoff))
        else:
            graph.add((potential, RDF.type, ASO.InteratomicPotential))

        graph.add((potential, ASO.hasReference, Literal(potential_doi)))
        graph.add((method, ASO.hasInteratomicPotential, potential))

    if software_id is not None:
        software_agent = URIRef(software_id)
        graph.add((software_agent, RDF.type, PROV.SoftwareAgent))
        graph.add((activity, PROV.wasAssociatedWith, software_agent))
        if software is not None:
            graph.add((software_agent, RDFS.label, Literal(software)))


def annotate_dft(graph,
                structure,
                id = None,
                pressure = None,
                software_id = None,
                software = None,
                ):
    """
    Annotate a given structure with MD simulation details
    """
    main_id = str(uuid.uuid4())
    if id is not None:
        main_id = id
        
    sample = structure.sample
    graph.add((sample, RDF.type, PROV.Entity))

    activity = URIRef(f'activity:{main_id}')
    graph.add((activity, RDF.type, PROV.Activity))
    graph.add((activity, RDF.type, ASO.StructureOptimization))

    method = URIRef(f'method:{main_id}')
    graph.add((method, RDF.type, ASO.DensityFunctionalTheory))
    graph.add((activity, ASO.hasMethod, method))
    graph.add((activity, ASO.hasRelaxationDOF, ASO.AtomicPosition))

    if pressure is None:
        pass
    elif np.isscalar(pressure):
        graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
    else: 
        #check if pressure is hydrostatic or not
        axial_all_alike = None not in pressure[:3] and np.allclose(
            pressure[:3], pressure[0]
        )
        shear_all_none = all(p is None for p in pressure[3:])
        shear_all_zero = None not in pressure[3:] and np.allclose(pressure[3:], 0)
        hydrostatic = axial_all_alike and (shear_all_none or shear_all_zero)
        if hydrostatic:
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
        else:
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
            graph.add((activity, ASO.hasRelaxationDOF, ASO.CellShape))


    graph.add((sample, PROV.wasGeneratedBy, activity))

    if software_id is not None:
        software_agent = URIRef(software_id)
        graph.add((software_agent, RDF.type, PROV.SoftwareAgent))
        graph.add((activity, PROV.wasAssociatedWith, software_agent))
        if software is not None:
            graph.add((software_agent, RDFS.label, Literal(software)))


def _get_inherited_properties(kg, from_sample, to_sample):
    #Here we need to add inherited info: CalculatedProperties will be lost
    #Defects will be inherited
    #add vac stuff
    material = list([k[2] for k in kg.graph.triples((from_sample, CMSO.hasMaterial, None))])[0]
    defects = list([x[2] for x in kg.graph.triples((material, CMSO.hasDefect, None))])
    #now for each defect we copy add this to the final sample
    final_material = list([k[2] for k in kg.graph.triples((to_sample, CMSO.hasMaterial, None))])[0]
    
    for defect in defects:
        new_defect = URIRef(defect.toPython())
        kg.graph.add((final_material, CMSO.hasDefect, new_defect))
        #now fetch all defect based info
        for triple in kg.graph.triples((defect, None, None)):
            kg.graph.add((new_defect, triple[1], triple[2]))
    
    #now add the special props for vacancy
    initial_simcell = kg.graph.value(from_sample, CMSO.hasSimulationCell)
    final_simcell = kg.graph.value(to_sample, CMSO.hasSimulationCell) 
    for triple in kg.graph.triples((initial_simcell, PODO.hasVacancyConcentration, None)):
        kg.graph.add((final_simcell, triple[1], triple[2]))
    for triple in kg.graph.triples((initial_simcell, PODO.hasNumberOfVacancies, None)):
        kg.graph.add((final_simcell, triple[1], triple[2]))

def _get_lattice_properties(kg, from_sample, to_sample):
    material = list([k[2] for k in kg.graph.triples((from_sample, CMSO.hasMaterial, None))])[0]
    crystal_structure = kg.graph.value(material, CMSO.hasStructure)
    altname = kg.graph.value(crystal_structure, CMSO.hasAltName)

    #add this to new structure
    final_material = list([k[2] for k in kg.graph.triples((to_sample, CMSO.hasMaterial, None))])[0]
    final_crystal_structure = kg.graph.value(final_material, CMSO.hasStructure)
    kg.add((final_crystal_structure, CMSO.hasAltName, altname))

    #space group
    space_group = kg.graph.value(crystal_structure, CMSO.hasSpaceGroup)
    final_space_group = kg.graph.value(final_crystal_structure, CMSO.hasSpaceGroup)
    for triple in kg.graph.triples((space_group, None, None)):
        kg.graph.add((final_space_group, triple[1], triple[2]))

    #unit cell
    unit_cell = kg.graph.value(crystal_structure, CMSO.hasUnitCell)
    bv = kg.graph.value(unit_cell, CMSO.hasBravaisLattice)

    final_unit_cell = kg.graph.value(final_crystal_structure, CMSO.hasUnitCell)
    kg.graph.add((final_unit_cell, CMSO.hasBravaisLattice, bv))

    #lattice parameter
    lattice_parameter = kg.graph.value(unit_cell, CMSO.hasLatticeParameter)
    final_lattice_parameter = kg.graph.value(final_unit_cell, CMSO.hasLatticeParameter)
    for triple in kg.graph.triples((lattice_parameter, None, None)):
        kg.graph.add((final_lattice_parameter, triple[1], triple[2]))
    
    #lattice angle
    lattice_angle = kg.graph.value(unit_cell, CMSO.hasAngle)
    final_lattice_angle = kg.graph.value(final_unit_cell, CMSO.hasAngle)
    for triple in kg.graph.triples((lattice_angle, None, None)):
        kg.graph.add((final_lattice_angle, triple[1], triple[2]))


def add_derived_structure(kg, initial_sample, final_sample):
    _get_inherited_properties(kg, initial_sample, final_sample)
    _get_lattice_properties(kg, initial_sample, final_sample)

    kg.add((initial_sample, RDF.type, PROV.Entity))
    kg.add((final_sample, RDF.type, PROV.Entity))
    kg.add((final_sample, PROV.wasDerivedFrom, initial_sample))

def add_method(kg, mdict):
    main_id = mdict['id']
    activity = URIRef(f'activity:{main_id}')
    kg.add((activity, RDF.type, PROV.Activity))

    if len(mdict['dof']) == 0:
        kg.add((activity, RDF.type, ASO.RigidEnergyCalculation))
    else:
        kg.add((activity, RDF.type, ASO.StructureOptimization))

    method = URIRef(f'method:{main_id}')
    if mdict['method'] == 'MolecularStatics':
        kg.add((method, RDF.type, ASO.MolecularStatics))
    elif mdict['method'] == 'MolecularDynamics':
        kg.add((method, RDF.type, ASO.MolecularDynamics))
    kg.add((activity, ASO.hasMethod, method))

    for dof in mdict['dof']:
        kg.add((activity, ASO.hasRelaxationDOF, getattr(ASO, dof)))

    kg.add((method, ASO.hasStatisticalEnsemble, getattr(ASO, mdict['ensemble'])))


    #add temperature if needed
    if mdict['temperature'] is not None:
        temperature = URIRef(f'temperature:{main_id}')
        kg.add((temperature, RDF.type, ASO.InputParameter))
        kg.add((temperature, RDFS.label, Literal('temperature', datatype=XSD.string)))
        kg.add((activity, ASO.hasInputParameter, temperature))
        kg.add((temperature, ASO.hasValue, Literal(mdict['temperature'], datatype=XSD.float)))
        kg.add((temperature, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))

    if mdict['pressure'] is not None:
        pressure = URIRef(f'pressure:{main_id}')
        kg.add((pressure, RDF.type, ASO.InputParameter))
        kg.add((pressure, RDFS.label, Literal('pressure', datatype=XSD.string)))
        kg.add((activity, ASO.hasInputParameter, pressure))
        kg.add((pressure, ASO.hasValue, Literal(mdict['pressure'], datatype=XSD.float)))
        kg.add((pressure, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/GigaPA')))

