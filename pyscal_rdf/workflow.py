"""
Workflows aspects for non-automated annotation of structures.

"""

from pyscal_rdf.rdfsystem import System
from pyscal_rdf import StructureGraph
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

