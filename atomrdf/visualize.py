import graphviz
import logging
import os
import re
import xml.etree.ElementTree as ET
from rdflib import BNode, URIRef, Namespace, Literal, RDF
from rdflib import Literal as RDFLiteral
import uuid
import json

logger = logging.getLogger(__name__)


# ─── GEXF export helpers ──────────────────────────────────────────────────────

#: RGB colours for each semantic category used in the GEXF / Gephi export.
GEXF_CATEGORY_COLORS = {
    "Sample":      (224, 123,  57),   # warm orange
    "Material":    (155,  89, 182),   # purple
    "Structure":   ( 41, 128, 185),   # blue
    "Element":     ( 39, 174,  96),   # green
    "Calculation": (192,  57,  43),   # red
    "Potential":   (243, 156,  18),   # gold
    "Property":    ( 22, 160, 133),   # teal
    "Literal":     (189, 195, 199),   # light grey
    "Other":       (149, 165, 166),   # grey
}

# exact rdf:type local-name → semantic category
_TYPE_CATEGORY_EXACT = {
    # Sample
    "AtomicScaleSample":              "Sample",
    # Material
    "CrystallineMaterial":            "Material",
    "AmorphousStructure":             "Material",
    "FluidMaterial":                  "Material",
    "CrystallineFluid":               "Material",
    # Structure — crystal geometry
    "CrystalStructure":               "Structure",
    "UnitCell":                       "Structure",
    "CrystalUnitCell":                "Structure",
    "BravaisLattice":                 "Structure",
    "SimulationCell":                 "Structure",
    "SimulationCellAngle":            "Structure",
    "SimulationCellLength":           "Structure",
    "SimulationCellVector":           "Structure",
    "LatticeParameter":               "Structure",
    "LatticeAngle":                   "Structure",
    # Structure — defects
    "GrainBoundary":                  "Structure",
    "TiltGrainBoundary":              "Structure",
    "TwistGrainBoundary":             "Structure",
    "SymmetricalTiltGrainBoundary":   "Structure",
    "Vacancy":                        "Structure",
    "InterstitialImpurity":           "Structure",
    "SubstitutionalImpurity":         "Structure",
    "DefectComplex":                  "Structure",
    # Element
    "ChemicalElement":                "Element",
    "Species":                        "Element",
    "ChemicalSpecies":                "Element",
    # Calculation — simulation types
    "MolecularStatics":               "Calculation",
    "MolecularDynamics":              "Calculation",
    "DensityFunctionalTheory":        "Calculation",
    "Simulation":                     "Calculation",
    "EnergyCalculation":              "Calculation",
    "GeneralizedGradientApproximation": "Calculation",
    "LocalDensityApproximation":      "Calculation",
    "InputParameter":                 "Calculation",
    "SoftwareAgent":                  "Calculation",
    # Calculation — relaxation types & ensembles
    "AtomicPositionRelaxation":       "Calculation",
    "CellShapeRelaxation":            "Calculation",
    "CellVolumeRelaxation":           "Calculation",
    "MicrocanonicalEnsemble":         "Calculation",
    "IsothermalIsobaricEnsemble":     "Calculation",
    "CanonicalEnsemble":              "Calculation",
    "GrandCanonicalEnsemble":         "Calculation",
    # Calculation — operations
    "DeleteAtom":                     "Calculation",
    "SubstituteAtom":                 "Calculation",
    "AddAtom":                        "Calculation",
    "Multiplication":                 "Calculation",
    "Division":                       "Calculation",
    "Subtraction":                    "Calculation",
    "Addition":                       "Calculation",
    # Potential
    "InteratomicPotential":           "Potential",
    "EmbeddedAtomModel":              "Potential",
    "ModifiedEmbeddedAtomModel":      "Potential",
    "PairPotential":                  "Potential",
    "FinnisSinclairPotential":        "Potential",
    "LennardJonesPotential":          "Potential",
    "MachineLearningPotential":       "Potential",
    "MACE":                           "Potential",
    "ACE":                            "Potential",
    "GRACE":                          "Potential",
    "NNPotential":                    "Potential",
    "DeePMD":                         "Potential",
    # Property
    "CalculatedProperty":             "Property",
    "GrainBoundaryEnergy":            "Property",
    "FormationEnergy":                "Property",
    "SurfaceEnergy":                  "Property",
    "VacancyFormationEnergy":         "Property",
    "SegregationEnergy":              "Property",
    "WorkOfSeparation":               "Property",
    "Volume":                         "Property",
    "TotalEnergy":                    "Property",
    "Energy":                         "Property",
    "Temperature":                    "Property",
    "Pressure":                       "Property",
    "Stress":                         "Property",
    "ElasticConstant":                "Property",
    "BulkModulus":                    "Property",
    "AtomAttribute":                  "Property",
}

# Regex that matches the UUID portion of atomRDF instance URIs such as
# ``sample:3f8a7b2c-1234-5678-abcd-ef0123456789_CrystalStructure``
_ATOMRDF_UUID_RE = re.compile(
    r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}_',
    re.IGNORECASE,
)


def _gexf_local_name(uri_str):
    """Return the local name of a URI (fragment after the last '#' or '/')."""
    for sep in ("#", "/"):
        if sep in uri_str:
            return uri_str.rsplit(sep, 1)[-1]
    return uri_str


def _gexf_term_key(term):
    """Return a hashable key that uniquely identifies an RDF term."""
    if isinstance(term, RDFLiteral):
        return f"lit\x00{term}\x00{term.datatype}"
    if isinstance(term, BNode):
        return f"bnode\x00{term}"
    return str(term)


def _gexf_label(term):
    """Return a short human-readable label for any RDF term."""
    if isinstance(term, RDFLiteral):
        s = str(term)
        return (s[:60] + "\u2026") if len(s) > 60 else s
    if isinstance(term, BNode):
        return str(term)[:24]
    return _gexf_local_name(str(term))


def _classify_by_local_name(local):
    """
    Return a semantic category for a URI local name, or ``None`` if unknown.

    Used for both rdf:type classification and URI-suffix extraction.
    """
    if local in _TYPE_CATEGORY_EXACT:
        return _TYPE_CATEGORY_EXACT[local]
    if "Sample" in local:
        return "Sample"
    if any(k in local for k in ("Potential", "ForceField")):
        return "Potential"
    if any(k in local for k in ("Property", "Energy", "Pressure",
                                "Stress", "Elasticity",
                                "Force", "Segregation", "Separation",
                                "Modulus", "Constant")):
        return "Property"
    if any(k in local for k in ("Statics", "Dynamics", "Functional",
                                "Simulation", "Calculation", "Activity",
                                "Theory", "Method", "Approximation")):
        return "Calculation"
    if any(k in local for k in ("Element", "Species")):
        return "Element"
    if any(k in local for k in ("Structure", "Lattice", "UnitCell", "Cell",
                                "GrainBoundary", "Grain", "Boundary",
                                "Defect", "Vacancy", "Impurity",
                                "Interstitial", "Substitutional")):
        return "Structure"
    if any(k in local for k in ("Material", "Crystalline", "Amorphous")):
        return "Material"
    return None


def _gexf_classify(term, type_map):
    """
    Return the semantic category string for *term*.

    Priority
    --------
    1. Literals → ``"Literal"``.
    2. ``rdf:type`` lookup — exact local-name match, then keyword match.
    3. URI suffix extraction — for ``sample:UUID_TypeSuffix`` URIs whose
       type is a Wikidata QID or is otherwise not in the type map.
    4. URI-scheme shortcuts (``simulation:``, ``property:``, …) for nodes
       that carry no explicit type.
    5. Fallback → ``"Other"``.
    """
    if isinstance(term, RDFLiteral):
        return "Literal"

    uri_str = str(term)
    key = _gexf_term_key(term)

    # ── 1. rdf:type classification (primary) ──────────────────────────────────
    for type_uri in type_map.get(key, ()):
        local = _gexf_local_name(str(type_uri))
        cat = _classify_by_local_name(local)
        if cat is not None:
            return cat

    # ── 2. URI suffix extraction ───────────────────────────────────────────────
    # atomRDF encodes the node's class in the URI suffix, e.g.
    # sample:3f8a..._CrystalStructure
    m = _ATOMRDF_UUID_RE.search(uri_str)
    if m:
        suffix = uri_str[m.end():]          # e.g. "CrystalStructure" or "SimulationCellVector_3"
        suffix = suffix.split("_")[0]       # drop any trailing index component
        cat = _classify_by_local_name(suffix)
        if cat is not None:
            return cat

    # ── 3. Direct local-name of the URI itself ────────────────────────────────
    # Catches pure ontology class URIs such as cmso:SimulationCellVector,
    # cmso:ChemicalElement, asmo:GrainBoundaryEnergy, etc., which appear as
    # the *objects* of rdf:type triples and have no type of their own.
    local = _gexf_local_name(uri_str)
    cat = _classify_by_local_name(local)
    if cat is not None:
        return cat

    # ── 4. URI-scheme shortcuts (fallback for untyped instance nodes) ─────────
    if "sample:" in uri_str:
        return "Sample"
    if any(tok in uri_str for tok in ("simulation:", "activity:", "operation:")):
        return "Calculation"
    if "property:" in uri_str:
        return "Property"
    if "potential:" in uri_str:
        return "Potential"

    return "Other"


def to_gexf(g, output_file, include_literals=False, positions=None, sizes=None,
            top_n_labels=None, label_overrides=None, top_label_uris=None,
            injected_type_map=None):
    """
    Export an RDF graph to GEXF format for visualisation in Gephi.

    Nodes are coloured by semantic category:

    =========== ====================== =========================================
    Category    Colour                 Covers
    =========== ====================== =========================================
    Sample      orange  ``#E07B39``    ``cmso:AtomicScaleSample`` instances
    Material    purple  ``#9B59B6``    Material description nodes
    Structure   blue    ``#2980B9``    Crystal-structure / unit-cell nodes
    Element     green   ``#27AE60``    Chemical element / species nodes
    Calculation red     ``#C0392B``    Simulation / activity nodes
    Potential   gold    ``#F39C12``    Interatomic potential nodes
    Property    teal    ``#16A085``    Calculated-property nodes
    Literal     l.grey  ``#BDC3C7``    RDF literal values
    Other       grey    ``#95A5A6``    Ontology terms & everything else
    =========== ====================== =========================================

    The ``viz:color`` attribute written into the GEXF file is read natively by
    Gephi and drives the default node colour.  The ``category`` attribute is
    also stored as a node attribute so it can be used in Gephi's *Partition*
    panel for colour/size adjustments after import.

    Parameters
    ----------
    g : rdflib.Graph
        The graph to serialise (plain, named, or conjunctive graph).
    output_file : str
        Destination path for the ``.gexf`` file.
    include_literals : bool, optional
        Whether to create a node for every literal value.  Default is
        ``False`` (drops literal nodes and their edges), which produces a
        cleaner resource-only graph that is easier to explore in Gephi.
    positions : dict, optional
        Mapping ``{uri_string: (x, y)}`` of pre-computed layout coordinates.
        Written as ``viz:position`` elements so Gephi uses them directly.
        When ``None`` no position attributes are written.
    sizes : dict, optional
        Mapping ``{uri_string: float}`` of pre-computed node sizes.
        Written as ``viz:size`` elements.  When ``None`` Gephi uses its default.
    top_n_labels : int, optional
        When set, only the *top_n_labels* highest-degree nodes keep a visible
        label string; all other nodes are exported with an empty label.  This
        is useful when opening in Gephi with "Show node labels" enabled —
        only the most connected nodes will display text.  When ``None`` (the
        default) every node keeps its label.
    label_overrides : dict, optional
        Mapping ``{uri_string: display_label}`` of explicit label replacements.
        Applied after ``_gexf_label()`` so any URI can be given a clean short
        name (e.g. ``{"http://www.vasp.at": "VASP"}``).

    Returns
    -------
    output_file : str
        The path of the file that was written.
    """
    GEXF_NS = "http://gexf.net/1.3"
    VIZ_NS  = "http://gexf.net/1.3/viz"

    ET.register_namespace("",    GEXF_NS)
    ET.register_namespace("viz", VIZ_NS)

    # ── 1. Build rdf:type map (term_key → set of type URIs) ───────────────────
    # If the caller pre-built the type_map from a richer graph (e.g. one that
    # still contains rdf:type triples even though those were stripped from g),
    # use it directly so node classification still works correctly.
    if injected_type_map is not None:
        type_map = injected_type_map
    else:
        type_map = {}
        for s, p, o in g:
            if p == RDF.type and isinstance(o, URIRef):
                k = _gexf_term_key(s)
                type_map.setdefault(k, set()).add(o)

    # ── 2. Collect unique nodes and directed edges ─────────────────────────────
    term_to_idx = {}   # term_key → integer index
    node_rows   = []   # [(label, category, uri_str), …]
    edge_rows   = []   # [(src_idx, tgt_idx, predicate_local_name), …]

    def _register(term):
        k = _gexf_term_key(term)
        if k not in term_to_idx:
            idx = len(node_rows)
            term_to_idx[k] = idx
            cat = _gexf_classify(term, type_map)
            uri = str(term)
            lbl = (label_overrides.get(uri) or _gexf_label(term)) if label_overrides else _gexf_label(term)
            node_rows.append((lbl, cat, uri))
        return term_to_idx[k]

    for s, p, o in g:
        if isinstance(o, RDFLiteral) and not include_literals:
            continue
        s_idx = _register(s)
        o_idx = _register(o)
        edge_rows.append((s_idx, o_idx, _gexf_local_name(str(p))))

    # ── 2b. Blank labels for nodes not in the labelled set ────────────────────
    # top_label_uris takes priority: an explicit set of URIs to keep labelled.
    # Fallback: top_n_labels uses internal degree to pick the top-N.
    if top_label_uris is not None:
        node_rows = [
            (lbl if uri in top_label_uris else "", cat, uri)
            for lbl, cat, uri in node_rows
        ]
    elif top_n_labels is not None:
        degree = [0] * len(node_rows)
        for src, tgt, _ in edge_rows:
            degree[src] += 1
            degree[tgt] += 1
        threshold = sorted(degree, reverse=True)[min(top_n_labels, len(degree)) - 1]
        top_idx = set(i for i, d in enumerate(degree) if d >= threshold)
        if len(top_idx) > top_n_labels:
            ranked = sorted(top_idx, key=lambda i: degree[i], reverse=True)
            top_idx = set(ranked[:top_n_labels])
        node_rows = [
            (lbl if i in top_idx else "", cat, uri)
            for i, (lbl, cat, uri) in enumerate(node_rows)
        ]

    # ── 3. Build XML tree ──────────────────────────────────────────────────────
    root = ET.Element(f"{{{GEXF_NS}}}gexf", version="1.3")

    meta_el = ET.SubElement(root, f"{{{GEXF_NS}}}meta")
    ET.SubElement(meta_el, f"{{{GEXF_NS}}}creator").text = "atomRDF"
    ET.SubElement(meta_el, f"{{{GEXF_NS}}}description").text = (
        "Knowledge graph exported from atomRDF"
    )

    graph_el = ET.SubElement(
        root, f"{{{GEXF_NS}}}graph", mode="static", defaultedgetype="directed"
    )

    # Node attribute declarations
    node_attrs_el = ET.SubElement(
        graph_el, f"{{{GEXF_NS}}}attributes", {"class": "node"}
    )
    ET.SubElement(node_attrs_el, f"{{{GEXF_NS}}}attribute",
                  id="cat", title="category", type="string")
    ET.SubElement(node_attrs_el, f"{{{GEXF_NS}}}attribute",
                  id="uri", title="uri", type="string")

    # Edge attribute declarations
    edge_attrs_el = ET.SubElement(
        graph_el, f"{{{GEXF_NS}}}attributes", {"class": "edge"}
    )
    ET.SubElement(edge_attrs_el, f"{{{GEXF_NS}}}attribute",
                  id="pred", title="predicate", type="string")

    # Nodes
    nodes_el = ET.SubElement(graph_el, f"{{{GEXF_NS}}}nodes")
    for idx, (lbl, cat, uri) in enumerate(node_rows):
        r, gv, b = GEXF_CATEGORY_COLORS.get(cat, (149, 165, 166))
        node_el = ET.SubElement(
            nodes_el, f"{{{GEXF_NS}}}node", id=str(idx), label=lbl
        )
        attvals_el = ET.SubElement(node_el, f"{{{GEXF_NS}}}attvalues")
        ET.SubElement(attvals_el, f"{{{GEXF_NS}}}attvalue",
                      {"for": "cat", "value": cat})
        ET.SubElement(attvals_el, f"{{{GEXF_NS}}}attvalue",
                      {"for": "uri", "value": uri})
        ET.SubElement(node_el, f"{{{VIZ_NS}}}color",
                      r=str(r), g=str(gv), b=str(b), a="1.0")
        if positions is not None and uri in positions:
            px, py = positions[uri]
            ET.SubElement(node_el, f"{{{VIZ_NS}}}position",
                          x=f"{float(px):.4f}", y=f"{float(py):.4f}", z="0.0")
        if sizes is not None and uri in sizes:
            ET.SubElement(node_el, f"{{{VIZ_NS}}}size",
                          value=f"{float(sizes[uri]):.4f}")

    # Edges
    edges_el = ET.SubElement(graph_el, f"{{{GEXF_NS}}}edges")
    for edge_idx, (src, tgt, pred) in enumerate(edge_rows):
        edge_el = ET.SubElement(
            edges_el, f"{{{GEXF_NS}}}edge",
            id=str(edge_idx), source=str(src), target=str(tgt), label=pred,
        )
        attvals_el = ET.SubElement(edge_el, f"{{{GEXF_NS}}}attvalues")
        ET.SubElement(attvals_el, f"{{{GEXF_NS}}}attvalue",
                      {"for": "pred", "value": pred})

    # ── 4. Write file ──────────────────────────────────────────────────────────
    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(ET.tostring(root, encoding="unicode"))

    # ── 5. Summary ────────────────────────────────────────────────────────────
    cat_counts = {}
    for _, cat, _ in node_rows:
        cat_counts[cat] = cat_counts.get(cat, 0) + 1

    logger.info("Exported to '%s'", output_file)
    logger.info("  Nodes : %s", f"{len(node_rows):,}")
    logger.info("  Edges : %s", f"{len(edge_rows):,}")
    logger.info("  Node categories:")
    for cat in ("Sample", "Material", "Structure", "Element",
                "Calculation", "Potential", "Property", "Literal", "Other"):
        count = cat_counts.get(cat, 0)
        if count:
            rv, gv, bv = GEXF_CATEGORY_COLORS[cat]
            logger.info("    %-12s %6s  #%02X%02X%02X", cat, f"{count:,}", rv, gv, bv)

    return output_file


def get_title_from_BNode(x):
    return x.toPython()


def get_string_from_URI(x):
    """
    Extract a presentable string from URI.

    Parameters
    ----------
    x : rdflib.term.URIRef
        The URI object to extract the string from.

    Returns
    -------
    tuple
        A tuple containing the presentable string representation of the URI and its type.
        The string representation is the last part of the URI after splitting by '#' or '/'.
        The type can be either "URIRef" or "BNode".
    """
    raw = x.toPython()
    # first try splitting by #
    rawsplit = raw.split("#")
    if len(rawsplit) > 1:
        return rawsplit[-1], "URIRef"

    # try splitting by = for chebi values
    if "CHEBI" in raw:
        rawsplit = raw.split("=")
        rawsplit = rawsplit[-1].split(":")
        if len(rawsplit) > 1:
            return ".".join(rawsplit[-2:]), "URIRef"

    if "sample:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "activity:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "simulation:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "operation:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    if "property:" in raw:
        rawsplit = raw.split(":")
        if len(rawsplit) > 1:
            return "_".join(rawsplit), "BNode"

    # just a normal url split now
    rawsplit = raw.split("/")
    if len(rawsplit) > 1:
        return ".".join(rawsplit[-2:]), "URIRef"

    # none of the conditions worked, which means it's a hex string
    return raw, "BNode"


def parse_object(x):
    """
    Parse the given object and return its title and type.

    Parameters
    ----------
    x : RDF term
        The RDF term to parse.

    Returns
    -------
    tuple
        A tuple containing the title of the object and its type.

    """
    if isinstance(x, BNode):
        return get_title_from_BNode(x), "BNode"
    elif isinstance(x, URIRef):
        return get_string_from_URI(x)
    elif isinstance(x, Literal):
        return str(x.title()), "Literal"


styledict = {
    "BNode": {"color": "#ffe6ff", "shape": "box", "style": "filled"},
    "URIRef": {"color": "#ffffcc", "shape": "box", "style": "filled"},
    "Literal": {"color": "#e6ffcc", "shape": "ellipse", "style": "filled"},
}


def _switch_box(box):
    if box == "box":
        return "rectangle"
    # remember that only boxes will be used, circles no!


def _fix_id(string1, istype1):
    if istype1 == "Literal":
        id1 = str(uuid.uuid4())
    else:
        id1 = string1
    return id1


def visualize_graph(
    g,
    styledict=styledict,
    rankdir="TB",
    hide_types=False,
    workflow_view=False,
    sample_view=False,
    size=None,
    layout="dot",
):
    """
    Visualizes a graph using Graphviz.

    Parameters
    ----------
    g : dict
        The graph to visualize.
    styledict : dict, optional
        A dictionary containing styles for different types of nodes and edges. Default is `styledict`.
    rankdir : str, optional
        The direction of the graph layout. Default is "TB" (top to bottom).
    hide_types : bool, optional
        Whether to hide nodes with the "type" attribute. Default is False.
    workflow_view : bool, optional
        Whether to enable the workflow view. Default is False.
    sample_view : bool, optional
        Whether to enable the sample view. Default is False.
    size : str, optional
        The size of the graph. Default is None.
    layout : str, optional
        The layout algorithm to use. Default is "dot".

    Returns
    -------
    dot : graphviz.Digraph
        The graph visualization.
    """
    dot = graphviz.Digraph()

    dot.attr(
        rankdir=rankdir,
        style="filled",
        size=size,
        layout=layout,
        overlap="false",
    )

    for k in g:
        string1, istype1 = parse_object(k[0])
        string2, istype2 = parse_object(k[2])
        string3, istype = parse_object(k[1])

        plot = True

        if workflow_view:
            # we collapse sample information
            # if cmso.connector is found, only use it is it is cmso.hasCalculated
            # all sub sample props, indicated by sample_x_jsjsj will be ignored.
            green_list = ["hasCalculatedProperty", "wasCalculatedBy", "hasValue"]
            ssplit = string3.split(".")
            if len(ssplit) == 2:
                if (ssplit[0] == "cmso") and (ssplit[1] not in green_list):
                    plot = False
            if string3 == "subClassOf":
                plot = False
            ssplit = string2.split(".")
            if string3 == "type":
                if (ssplit[0] == "cmso") and (ssplit[1] not in ["CalculatedProperty"]):
                    plot = False
                if (ssplit[0] == "cmso") and (ssplit[1] in ["AtomicScaleSample"]):
                    dot.node(
                        string1,
                        label=string1,
                        shape=styledict[istype1]["shape"],
                        style=styledict[istype1]["style"],
                        color=styledict[istype1]["color"],
                        fontsize=styledict[istype1]["fontsize"],
                        fontname=styledict[istype1]["fontname"],
                    )
                    plot = False
        
        elif sample_view:
            green_list = ['wasDerivedFrom', 'wasGeneratedBy']
            if string3 not in green_list:
                plot = False
            

        if hide_types and (string3 == "type"):
            plot = False

        if not plot:
            continue

        if istype1 == "Literal":
            id1 = str(uuid.uuid4())
        else:
            id1 = string1
        dot.node(
            id1,
            label=string1,
            shape=styledict[istype1]["shape"],
            style=styledict[istype1]["style"],
            color=styledict[istype1]["color"],
            fontsize=styledict[istype1]["fontsize"],
            fontname=styledict[istype1]["fontname"],
        )

        if istype2 == "Literal":
            id2 = str(uuid.uuid4())
        else:
            id2 = string2
        dot.node(
            id2,
            label=string2,
            shape=styledict[istype2]["shape"],
            style=styledict[istype2]["style"],
            color=styledict[istype2]["color"],
            fontsize=styledict[istype2]["fontsize"],
            fontname=styledict[istype2]["fontname"],
        )

        dot.edge(
            id1,
            id2,
            color=styledict["edgecolor"],
            label=string3,
            fontsize=styledict[istype2]["fontsize"],
            fontname=styledict[istype2]["fontname"],
        )

    return dot

def _id(item):
    return str(item).replace(':', '_')

def visualize_provenance(
    prov,
    rankdir="TB",
    size=None,
    layout="dot",
):
    dot = graphviz.Digraph()
    dot.attr(
        rankdir=rankdir,
        style="filled",
        size=size,
        layout=layout,
        overlap="false",
    )
    #add all nodes
    for key in prov.keys():
        nid = _id(key)
        #if "activity" in key:
        dot.node(nid, label=prov[key]['label'], 
                shape='box', 
                color="#C9DAF8", 
                style="filled",
                fontname='Helvetica',
                fontsize='8')
        #else:
        #    dot.node(nid, label=prov[key]['label'], 
        #            shape='parallelogram', 
        #            color="#C9DAF8", 
        #            style="filled",
        #            fontname='Helvetica',
        #            fontsize='8')
    #add all edges
    for key, val in prov.items():
        if 'inputs' in val.keys():
            if val['operation'] == 'input_parameter':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='input_param', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            if val['operation'] == 'output_parameter':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='output_param', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            elif val['operation'] == 'sample_for_activity':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='input_sample', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            elif val['operation'] == 'sample_output':
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), _id(key), label='output_sample', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
            else:
                operation_id = str(uuid.uuid4())
                operation = dot.node(operation_id, label=val['operation'], 
                                    color="#E6B8AF", 
                                    shape='box', 
                                    style='filled',
                                    fontname='Helvetica',
                                    fontsize='8')
                for subkey, subval in val['inputs'].items():
                    dot.edge(_id(subval), operation_id, label='input', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
                dot.edge(operation_id, _id(key), label='output', 
                        color="#263238",
                        fontname='Helvetica',
                        fontsize='8')
    return dot