# utility methods
from atomrdf.namespace import CMSO
from rdflib import URIRef
from collections import Counter
import numpy as np


def get_material(kg, sample_id):
    return kg.value(URIRef(sample_id), CMSO.hasMaterial)
