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

    cmso = OntologyNetwork(os.path.join(file_location, "data/cmso.owl")) 
    pldo = OntologyNetwork(os.path.join(file_location, "data/pldo.owl")) 
    podo = OntologyNetwork(os.path.join(file_location, "data/podo.owl")) 
    asmo = OntologyNetwork(os.path.join(file_location, "data/asmo.owl")) 
    ldo = OntologyNetwork(os.path.join(file_location, "data/ldo.owl")) 
    cdco = OntologyNetwork(os.path.join(file_location, "data/cdco.owl")) 

    # combine them
    combo = cmso + pldo + podo + asmo + ldo + cdco

    # add namespaces
    #combo.add_namespace("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    #combo.add_namespace("rdfs", "http://www.w3.org/2000/01/rdf-schema#")

    # add extra terms for quering
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

    # interontology paths
    #CMSO -> PODO VACANCY
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "pldo:PlanarDefect"))
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "podo:Vacancy"))
    combo.add_path(("cmso:ComputationalSample", "podo:hasVacancyConcentration", "float"))
    combo.add_path(("cmso:ComputationalSample", "podo:hasNumberOfVacancies", "int"))
    
    #CMSO -> PODO IMPURITY
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "podo:SubstitutionalImpurity"))
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "podo:InterstitialImpurity"))
    #combo.add_path(("podo:InterstitialImpurity", "rdfs:label", "string"))
    combo.add_path(("cmso:ComputationalSample", "podo:hasNumberOfImpurityAtoms", "int"))
    combo.add_path(("cmso:ComputationalSample", "podo:hasImpurityConcentration", "float"))

    #CMSO -> LDO DISL paths
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "ldo:LineDefect"))
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "ldo:Dislocation"))

    #CMSO -> PODO DISL paths
    combo.add_path(("cmso:CrystallineMaterial", "cdco:hasCrystallographicDefect", "pldo:GrainBoundary"))
    
    # return
    return combo
