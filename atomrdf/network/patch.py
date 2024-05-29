"""
These are patches specifically designed for atomrdf.

These may or may not be implemented in the ontology. As it is implemented; it can be removed
from the patches
"""

import os


def patch_terms(iri, rn):
    """
    Remove functions as patching is done
    """
    # Term: hasSymbol
    # Ontology: CMSO
    # Reason: Range is not specified in the owl file.
    # This prevents owlready2 from reading in this property correctly.
    if iri == "http://purls.helmholtz-metadaten.de/cmso/hasSymbol":
        rn = ["str"]
    # Term: hasValue
    # Ontology: CMSO
    # Reason: Range is Literal(); however here we use this for number values, hence we can fix this.
    # See fn: `add_calculated_property`
    elif iri == "http://purls.helmholtz-metadaten.de/cmso/hasValue":
        rn = ["float", "double", "int", "str"]

    elif iri == "http://purls.helmholtz-metadaten.de/asmo/hasValue":
        rn = ["float", "double", "int", "str"]

    #elif iri == "http://purls.helmholtz-metadaten.de/asmo/hasUnit":
    #    rn = ["str"]

    #elif iri == "http://purls.helmholtz-metadaten.de/cmso/hasUnit":
    #    rn = ["str"]

    elif iri == "http://www.w3.org/2000/01/rdf-schema#label":
        rn = ["str"]

    elif iri == "http://purls.helmholtz-metadaten.de/cmso/hasReference":
        rn = ["str"]

    elif iri == 'http://purls.helmholtz-metadaten.de/cmso/hasSpaceGroupSymbol':
        rn = ["str"]

    return rn
