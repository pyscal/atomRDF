# utility methods
from atomrdf.namespace import CMSO
from rdflib import URIRef
from collections import Counter
import numpy as np


def get_material(kg, sample_id):
    if isinstance(sample_id, str):
        sample_id = URIRef(sample_id)
    return kg.value(sample_id, CMSO.hasMaterial)


def get_sample_id(sample_id):
    if isinstance(sample_id, str):
        return sample_id
    return str(sample_id.toPython())
