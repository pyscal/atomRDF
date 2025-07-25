# utility methods
from atomrdf.namespace import CMSO
from rdflib import URIRef


def get_material(kg, sample_id):
    return kg.value(URIRef(sample_id), CMSO.hasMaterial)
