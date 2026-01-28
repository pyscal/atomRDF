# utility methods
from atomrdf.namespace import CMSO
from rdflib import URIRef
from collections import Counter
import numpy as np


def get_simulation(kg, simulation_id):
    if isinstance(simulation_id, str):
        simulation_id = URIRef(simulation_id)
    return simulation_id


def get_material(kg, sample_id):
    if isinstance(sample_id, str):
        sample_id = URIRef(sample_id)
    return kg.value(sample_id, CMSO.hasMaterial)


def get_sample_id(sample_id):
    if isinstance(sample_id, str):
        return sample_id
    return str(sample_id.toPython())


def get_sample_object(sample_id):
    if isinstance(sample_id, str):
        sample_id = URIRef(sample_id)
    return sample_id


def toPython(item):
    if item is None:
        return None
    try:
        return item.toPython()
    except AttributeError:
        return item
