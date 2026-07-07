"""Ontology loading for atomRDF's ontology-based ("auto") query feature.

Builds the combined OCDO ontology network consumed by
:meth:`atomrdf.KnowledgeGraph.query` for ontology-term queries. This depends on
``tools4rdf``, distributed as the optional ``autoquery`` extra
(``pip install atomrdf[autoquery]``). It is imported lazily so the rest of
atomRDF works without ``tools4rdf`` installed.
"""


def read_ontology():
    try:
        from tools4rdf.network.parser import parse_ontology
        from tools4rdf.network.network import OntologyNetworkBase
    except ImportError as exc:
        raise ImportError(
            "Ontology-based querying needs the 'tools4rdf' package, which ships "
            "as the optional 'autoquery' extra. Install it with:\n"
            "    pip install atomrdf[autoquery]"
        ) from exc

    cmso = parse_ontology("https://purls.helmholtz-metadaten.de/cmso/")
    pldo = parse_ontology("https://purls.helmholtz-metadaten.de/cdos/pldo/")
    podo = parse_ontology("https://purls.helmholtz-metadaten.de/cdos/podo/")
    asmo = parse_ontology("https://purls.helmholtz-metadaten.de/asmo/")
    ldo = parse_ontology("https://purls.helmholtz-metadaten.de/cdos/ldo/")
    cdco = parse_ontology("https://purls.helmholtz-metadaten.de/cdos/cdco/")

    # now sum them up
    combo = cmso + cdco + pldo + podo + asmo + ldo
    combo.attributes["data_property"]["cmso:hasSymbol"].range.append("str")
    combo.attributes["data_property"]["asmo:hasValue"].range.extend(
        ["float", "double", "int", "str"]
    )

    # now combine the ontologies
    combo = OntologyNetworkBase(combo)

    # add sring labels as needed
    combo.add_namespace("rdfs", "http://www.w3.org/2000/01/rdf-schema#")

    combo.add_term(
        "http://www.w3.org/2000/01/rdf-schema#label",
        "data_property",
        delimiter="#",
        namespace="rdfs",
        rn=["str"],
    )

    combo.add_path(("asmo:CalculatedProperty", "rdfs:label", "string"))
    combo.add_path(("asmo:InputParameter", "rdfs:label", "string"))
    combo.add_path(("prov:SoftwareAgent", "rdfs:label", "string"))
    combo.add_path(("asmo:InteratomicPotential", "rdfs:label", "string"))

    return combo
