from tools4rdf.network.parser import OntoParser, parse_ontology
from tools4rdf.network.network import OntologyNetworkBase

import os

def read_ontology():
    file_location = os.path.dirname(__file__)

    cmso = parse_ontology(os.path.join(file_location, "data/cmso.owl")) 
    pldo = parse_ontology(os.path.join(file_location, "data/pldo.owl")) 
    podo = parse_ontology(os.path.join(file_location, "data/podo.owl")) 
    asmo = parse_ontology(os.path.join(file_location, "data/asmo.owl")) 
    ldo = parse_ontology(os.path.join(file_location, "data/ldo.owl")) 
    cdco = parse_ontology(os.path.join(file_location, "data/cdco.owl"))


    #now sum them up
    combo = cmso + cdco + pldo + podo + asmo + ldo
    combo.attributes['data_property']['cmso:hasSymbol'].range.append("str")
    combo.attributes['data_property']['asmo:hasValue'].range.extend(["float", "double", "int", "str"])

    #now combine the ontologies
    combo = OntologyNetworkBase(combo)

    #add sring labels as needed
    combo.add_namespace("rdfs", "http://www.w3.org/2000/01/rdf-schema#")

    combo.add_term(
        "http://www.w3.org/2000/01/rdf-schema#label",
        "data_property",
        delimiter="#",
        namespace="rdfs",
        rn = ['str']
    )

    combo.add_path(("asmo:CalculatedProperty", "rdfs:label", "string"))
    combo.add_path(("asmo:InputParameter", "rdfs:label", "string"))
    combo.add_path(("prov:SoftwareAgent", "rdfs:label", "string"))
    combo.add_path(("asmo:InteratomicPotential", "rdfs:label", "string"))

    return combo