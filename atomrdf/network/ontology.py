"""
Documentation for the ontology module.

Updates needed
--------------
ASMO module

- Add math operations, see sample.py for the complete list of items - at the moment, 4 cardinal operations,
but maybe more will be needed.

- StructureOperation: 
    - RotationOperation and associated vectors, see structure.py
    - TranslationOperation - no triples yet, see structure.py
    - ADD SHEAR OPERATION
- PointDefectCreation
    - DeleteAtom - to create vacancy - Now Just added as label to PROV activities
    - AddAtom - to create impurity
    - SubstituteAtom - to create impurity
"""

import os
from atomrdf.network.network import OntologyNetwork


def read_ontology():
    """
    Read in ontologies and perform necessary operations.

    Returns
    -------
    combo: OntologyNetwork, Combined ontology network.
    """
    # read in ontologies
    file_location = os.path.dirname(__file__).split("/")
    file_location = "/".join(file_location[:-1])

    cmso = OntologyNetwork(os.path.join(file_location, "data/cmso.owl")) #eb8a58a60af1759afc4115ae395ad62ddffb6844
    pldo = OntologyNetwork(os.path.join(file_location, "data/pldo.owl")) #d15d27712e3f64b405d75c70ad970c9d54ff0b51
    podo = OntologyNetwork(os.path.join(file_location, "data/podo.owl")) #6a74d511c5b78042e1cb7a6e76e948fa56de598e
    asmo = OntologyNetwork(os.path.join(file_location, "data/asmo.owl")) #c7e2da99b9126844f19f225c6a10cdb01aeb55e6
    ldo = OntologyNetwork(os.path.join(file_location, "data/ldo.owl")) #e23fa9930351787e701347878a3e1a0e3924d084

    # combine them
    combo = cmso + pldo + podo + asmo + ldo

    # add namespaces
    combo.add_namespace("prov", "http://www.w3.org/ns/prov#")
    combo.add_namespace("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    combo.add_namespace("rdfs", "http://www.w3.org/2000/01/rdf-schema#")

    # add extra terms for quering
    combo.add_term("http://www.w3.org/ns/prov#Entity", "class", delimiter="#")
    combo.add_term("http://www.w3.org/ns/prov#Activity", "class", delimiter="#")
    combo.add_term("http://www.w3.org/ns/prov#SoftwareAgent", "class", delimiter="#")
    combo.add_term(
        "http://www.w3.org/ns/prov#wasDerivedFrom", "object_property", delimiter="#"
    )
    combo.add_term(
        "http://www.w3.org/ns/prov#wasGeneratedBy", "object_property", delimiter="#"
    )
    combo.add_term(
        "http://www.w3.org/ns/prov#wasAssociatedWith", "object_property", delimiter="#"
    )
    combo.add_term(
        "http://www.w3.org/ns/prov#actedOnBehalfOf", "object_property", delimiter="#"
    )
    combo.add_term(
        "http://www.w3.org/2000/01/rdf-schema#label",
        "data_property",
        delimiter="#",
        namespace="rdfs",
        rn = ['str']
    )
    combo.add_term(
        "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
        "object_property",
        delimiter="#",
        namespace="rdf",
    )

    # add paths

    # General fixes
    combo.add_path(("cmso:CrystalStructure", "cmso:hasAltName", "string"))

    # interontology paths
    #CMSO -> PODO VACANCY
    combo.add_path(("cmso:Material", "cmso:hasDefect", "pldo:PlanarDefect"))
    combo.add_path(("cmso:Material", "cmso:hasDefect", "podo:Vacancy"))
    combo.add_path(("cmso:AtomicScaleSample", "podo:hasVacancyConcentration", "float"))
    combo.add_path(("cmso:AtomicScaleSample", "podo:hasNumberOfVacancies", "int"))
    
    #CMSO -> PODO IMPURITY
    combo.add_path(("cmso:Material", "cmso:hasDefect", "podo:SubstitutionalImpurity"))
    combo.add_path(("cmso:Material", "cmso:hasDefect", "podo:InterstitialImpurity"))
    combo.add_path(("podo:InterstitialImpurity", "rdfs:label", "string"))
    combo.add_path(("cmso:AtomicScaleSample", "podo:hasNumberOfImpurityAtoms", "int"))
    combo.add_path(("cmso:AtomicScaleSample", "podo:hasImpurityConcentration", "float"))

    #CMSO -> LDO DISL paths
    combo.add_path(("cmso:Material", "cmso:hasDefect", "ldo:LineDefect"))
    combo.add_path(("cmso:Material", "cmso:hasDefect", "ldo:Dislocation"))


    combo.add_path(
        ("cmso:ComputationalSample", "prov:wasDerivedFrom", "cmso:ComputationalSample")
    )
    combo.add_path(("cmso:ComputationalSample", "rdf:type", "prov:Entity"))
    combo.add_path(("cmso:AtomicScaleSample", "prov:wasGeneratedBy", "prov:Activity"))
    combo.add_path(("asmo:EnergyCalculation", "rdf:type", "prov:Activity"))
    combo.add_path(
        ("asmo:EnergyCalculation", "prov:wasAssociatedWith", "prov:SoftwareAgent")
    )
    combo.add_path(
        (
            "cmso:ComputationalSample",
            "prov:wasGeneratedBy",
            "asmo:EnergyCalculation",
        )
    )
    
    combo.add_path(("cmso:CalculatedProperty", "asmo:hasValue", "float"))
    #combo.add_path(("cmso:CalculatedProperty", "asmo:hasUnit", "string"))
    combo.add_path(("cmso:CalculatedProperty", "rdfs:label", "string"))
    #how to handle units?

    #input parameters for ASMO
    combo.add_path(("asmo:InputParameter", "rdfs:label", "string"))
    #combo.add_path(("cmso:CalculatedProperty", "asmo:hasUnit", "string"))
    combo.add_path(("asmo:InputParameter", "asmo:hasValue", "float"))

    #software agent
    combo.add_path(("prov:SoftwareAgent", "rdfs:label", "string"))
    combo.add_path(("asmo:InteratomicPotential", "cmso:hasReference", "string"))
    combo.add_path(("asmo:InteratomicPotential", "rdfs:label", "string"))

    # return
    return combo
