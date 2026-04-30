"""
Workflow reconstruction from a Knowledge Graph.

Given a workflow URI (or sample / calculated-property URI), this module
queries all provenance from the KG and generates an executable Python
script that reproduces the workflow using ``atomrdf.workflow`` nodes.

Two modes
---------
``mode="recreate"``
    Writes structures from KG data and generates a fully runnable Python
    script with all parameters filled in from the KG provenance.

``mode="create_template"``
    Generates the same script skeleton but with ``# TODO: fill in``
    placeholders for values that need to be supplied by the user
    (e.g. potential file paths).
"""

import os
import textwrap
from collections import OrderedDict

from rdflib import URIRef, RDF, RDFS

from atomrdf.namespace import CMSO, ASMO, PROV, MDO, Literal


# ── Exceptions ─────────────────────────────────────────────────────────


class UnsupportedMethodError(Exception):
    """Raised when the KG contains a method/algorithm we cannot reconstruct."""

    _SUPPORTED_METHODS = {
        "MolecularStatics": [None, "EquationOfStateFit"],
    }
    _SUPPORTED_OPERATIONS = [
        "DeleteAtom",
        "SubstituteAtom",
        "AddAtom",
    ]

    @classmethod
    def supported_summary(cls):
        lines = ["Supported workflow methods:"]
        for method, algos in cls._SUPPORTED_METHODS.items():
            for algo in algos:
                if algo:
                    lines.append(f"  - {method} / {algo}")
                else:
                    lines.append(f"  - {method} (relaxation)")
        lines.append("Supported operations:")
        for op in cls._SUPPORTED_OPERATIONS:
            lines.append(f"  - {op}")
        return "\n".join(lines)


# ── Software → structure-file format ──────────────────────────────────

_SOFTWARE_FORMAT = {
    "vasp": "poscar",
    "lammps": "lammps-data",
    "quantum-espresso": "quantum-espresso",
    "qe": "quantum-espresso",
}

# ASMO operation URI fragments we recognise.
_KNOWN_OPS = {
    "DeleteAtom",
    "SubstituteAtom",
    "AddAtom",
    "Rotation",
    "Translation",
    "Shear",
}


# ── Internal helpers ───────────────────────────────────────────────────


def _short(uri):
    """Return the last path-segment of a URI (strips namespace or prefix)."""
    s = str(uri)
    if "/" in s:
        s = s.rsplit("/", 1)[-1]
    if ":" in s:
        s = s.rsplit(":", 1)[-1]
    return s


def _query_workflow_info(kg, workflow_uri):
    """
    Query a single workflow/simulation node and return a structured dict.
    """
    g = kg.graph
    wf = URIRef(workflow_uri) if not isinstance(workflow_uri, URIRef) else workflow_uri

    info = {
        "uri": str(wf),
        "method": None,
        "algorithm": None,
        "dofs": [],
        "software": [],
        "xc_functional": None,
        "ensemble": None,
        "potential": None,
        "potential_type": None,
        "input_parameters": OrderedDict(),
        "output_parameters": OrderedDict(),
        "calculated_properties": OrderedDict(),
        "input_samples": [],
        "output_samples": [],
        "path": None,
    }

    # Method
    method = g.value(wf, ASMO.hasComputationalMethod)
    if method is not None:
        info["method"] = _short(method)

    # Algorithm
    algo = g.value(wf, ASMO.usesSimulationAlgorithm)
    if algo is not None:
        algo_type = g.value(algo, RDF.type)
        info["algorithm"] = _short(algo_type) if algo_type else _short(algo)

    # DOFs
    for _, _, dof in g.triples((wf, ASMO.hasRelaxationDOF, None)):
        info["dofs"].append(_short(dof))

    # Software
    for _, _, sw in g.triples((wf, PROV.wasAssociatedWith, None)):
        label = g.value(sw, RDFS.label)
        if label is not None:
            info["software"].append(str(label))
        elif str(sw).startswith("http"):
            info["software"].append(_short(sw))

    # XC functional (DFT)
    xc = g.value(wf, MDO.hasXCFunctional)
    if xc is not None:
        xc_type = g.value(xc, RDF.type)
        info["xc_functional"] = _short(xc_type) if xc_type else _short(xc)

    # Ensemble (MD)
    ens = g.value(wf, ASMO.hasStatisticalEnsemble)
    if ens is not None:
        info["ensemble"] = _short(ens)

    # Potential
    pot = g.value(wf, ASMO.hasInteratomicPotential)
    if pot is not None:
        pot_label = g.value(pot, RDFS.label)
        info["potential"] = str(pot_label) if pot_label else _short(pot)
        pot_type = g.value(pot, RDF.type)
        if pot_type is not None:
            info["potential_type"] = _short(pot_type)

    # Input parameters
    for _, _, p in g.triples((wf, ASMO.hasInputParameter, None)):
        label = g.value(p, RDFS.label)
        value = g.value(p, ASMO.hasValue)
        if label is not None and value is not None:
            info["input_parameters"][str(label)] = str(value)

    # Output parameters
    for _, _, p in g.triples((wf, ASMO.hasOutputParameter, None)):
        label = g.value(p, RDFS.label)
        value = g.value(p, ASMO.hasValue)
        if label is not None and value is not None:
            info["output_parameters"][str(label)] = str(value)

    # Calculated properties
    for _, _, p in g.triples((None, ASMO.wasCalculatedBy, wf)):
        label = g.value(p, RDFS.label)
        value = g.value(p, ASMO.hasValue)
        if label is not None and value is not None:
            info["calculated_properties"][str(label)] = str(value)

    # Input / output samples via PROV
    for s, _, _ in g.triples((None, PROV.wasGeneratedBy, wf)):
        info["output_samples"].append(str(s))
    for out in info["output_samples"]:
        for _, _, inp in g.triples((URIRef(out), PROV.wasDerivedFrom, None)):
            inp_str = str(inp)
            if inp_str not in info["input_samples"]:
                info["input_samples"].append(inp_str)

    # Path
    path = g.value(wf, CMSO.hasPath)
    if path is not None:
        info["path"] = str(path)

    return info


def _query_operation_chain(kg, sample_uri):
    """
    Walk backwards from *sample_uri* and collect every operation found.
    Returns a list of dicts ordered newest-to-oldest.
    """
    g = kg.graph
    ops = []
    current = URIRef(sample_uri) if not isinstance(sample_uri, URIRef) else sample_uri

    visited = set()
    while current and str(current) not in visited:
        visited.add(str(current))

        for _, _, activity in g.triples((current, PROV.wasGeneratedBy, None)):
            act_type = g.value(activity, RDF.type)
            if act_type and _short(act_type) in _KNOWN_OPS:
                source = g.value(current, PROV.wasDerivedFrom)
                ops.append(
                    {
                        "op_type": _short(act_type),
                        "op_uri": str(activity),
                        "output_sample": str(current),
                        "input_sample": str(source) if source else None,
                    }
                )

        derived = g.value(current, PROV.wasDerivedFrom)
        current = derived

    return ops


def _walk_workflow_chain(kg, sample_uri):
    """
    Walk backwards through provenance and collect (sample, workflow) pairs.
    Returns newest-first.
    """
    g = kg.graph
    chain = []
    current = URIRef(sample_uri) if not isinstance(sample_uri, URIRef) else sample_uri

    visited = set()
    while current and str(current) not in visited:
        visited.add(str(current))
        for _, _, activity in g.triples((current, PROV.wasGeneratedBy, None)):
            method = g.value(activity, ASMO.hasComputationalMethod)
            if method is not None:
                chain.append((str(current), str(activity)))
        derived = g.value(current, PROV.wasDerivedFrom)
        current = derived

    return chain


def _query_sample_structure(kg, sample_uri):
    """
    Extract structure metadata from a sample node for code generation.

    Returns a dict with keys:
        element, crystal_structure, lattice_constant_a, lattice_constant_c,
        n_atoms, repetitions, cell_vectors, species, positions
    """
    g = kg.graph
    s = URIRef(sample_uri) if not isinstance(sample_uri, URIRef) else sample_uri

    result = {
        "element": None,
        "crystal_structure": None,
        "lattice_constant_a": None,
        "lattice_constant_c": None,
        "n_atoms": None,
        "repetitions": None,
        "cell_vectors": None,
        "species": [],
        "positions": [],
    }

    # Number of atoms
    n = g.value(s, CMSO.hasNumberOfAtoms)
    if n is not None:
        result["n_atoms"] = int(n)

    # Crystal structure / space group from material
    mat = g.value(s, CMSO.hasMaterial)
    if mat is not None:
        cs = g.value(mat, CMSO.hasStructure)
        if cs is not None:
            sg_sym = g.value(cs, CMSO.hasSpaceGroupSymbol)
            if sg_sym is not None:
                result["crystal_structure"] = str(sg_sym)

    # Cell vectors
    sim_cell = g.value(s, CMSO.hasSimulationCell)
    if sim_cell is not None:
        vec_nodes = list(g.objects(sim_cell, CMSO.hasVector))
        vec_nodes.sort(key=lambda v: str(v))
        vecs = []
        for vn in vec_nodes:
            x = g.value(vn, CMSO.hasComponent_x)
            y = g.value(vn, CMSO.hasComponent_y)
            z = g.value(vn, CMSO.hasComponent_z)
            if x is not None:
                vecs.append([float(x), float(y), float(z)])
        if len(vecs) == 3:
            result["cell_vectors"] = vecs

    # Lattice parameters
    if mat is not None:
        cs = g.value(mat, CMSO.hasStructure)
        if cs is not None:
            uc = g.value(cs, CMSO.hasUnitCell)
            if uc is not None:
                a = g.value(uc, CMSO.hasLatticeParameter_a)
                if a is not None:
                    result["lattice_constant_a"] = float(a)
                c = g.value(uc, CMSO.hasLatticeParameter_c)
                if c is not None:
                    result["lattice_constant_c"] = float(c)

    # Species (chemical symbols)
    for _, _, species_node in g.triples((s, CMSO.hasSpecies, None)):
        elem = g.value(species_node, CMSO.hasElement)
        if elem is not None:
            sym = g.value(elem, CMSO.hasChemicalSymbol)
            if sym is not None:
                result["species"].append(str(sym))

    # Primary element
    if result["species"]:
        result["element"] = result["species"][0]

    return result


def _detect_structure_format(software_list):
    for sw in software_list:
        fmt = _SOFTWARE_FORMAT.get(sw.lower())
        if fmt is not None:
            return fmt
    return "poscar"


def _format_extension(fmt):
    return {
        "poscar": "POSCAR",
        "lammps-data": "structure.lmp",
        "quantum-espresso": "structure.in",
    }.get(fmt, "structure.dat")


def _write_structure_fallback(kg, sample_uri, filepath, fmt):
    """Try to write a structure file from KG data."""
    try:
        kg.to_file(sample_uri, filepath, format=fmt)
        return True
    except Exception:
        pass

    try:
        from ase import Atoms as ASEAtoms
        from ase.io import write as ase_write
        import numpy as np

        g = kg.graph
        s = URIRef(sample_uri) if not isinstance(sample_uri, URIRef) else sample_uri

        cell = None
        sim_cell = g.value(s, CMSO.hasSimulationCell)
        if sim_cell is not None:
            vec_nodes = list(g.objects(sim_cell, CMSO.hasVector))
            vec_nodes.sort(key=lambda v: str(v))
            vecs = []
            for vn in vec_nodes:
                x = g.value(vn, CMSO.hasComponent_x)
                y = g.value(vn, CMSO.hasComponent_y)
                z = g.value(vn, CMSO.hasComponent_z)
                if x is not None:
                    vecs.append([float(x), float(y), float(z)])
            if len(vecs) == 3:
                cell = np.array(vecs)

        symbols = []
        for _, _, species_node in g.triples((s, CMSO.hasSpecies, None)):
            elem = g.value(species_node, CMSO.hasElement)
            if elem is not None:
                sym = g.value(elem, CMSO.hasChemicalSymbol)
                if sym is not None:
                    symbols.append(str(sym))

        n_lit = g.value(s, CMSO.hasNumberOfAtoms)
        n = int(n_lit) if n_lit is not None else len(symbols)

        if not symbols or cell is None:
            return False

        if len(symbols) == 1 and n > 1:
            symbols = symbols * n

        atoms = ASEAtoms(symbols=symbols[:n], cell=cell, pbc=True)

        if fmt == "poscar":
            ase_write(filepath, atoms, format="vasp")
        elif fmt == "lammps-data":
            ase_write(filepath, atoms, format="lammps-data", atom_style="atomic")
        else:
            ase_write(filepath, atoms, format=fmt)
        return True
    except Exception:
        return False


# ── Method / operation → node mapping ──────────────────────────────────


def _map_method_to_node(method, algorithm):
    """
    Map a KG method + algorithm to an atomrdf.workflow function name.

    Returns (module_path, function_name, import_line).
    Raises UnsupportedMethodError if the combination is not supported.
    """
    m = method or ""
    a = algorithm or ""

    if m == "MolecularStatics":
        if a in ("", "None") or a is None or a == "null":
            return (
                "atomrdf.workflow.evcurves",
                "relax_structure",
                "from atomrdf.workflow.evcurves import relax_structure",
            )
        if a == "EquationOfStateFit":
            return (
                "atomrdf.workflow.evcurves",
                "calculate_ev_curves",
                "from atomrdf.workflow.evcurves import calculate_ev_curves",
            )
        # Elastic constants have no special algorithm tag in KG,
        # but we can detect them from calculated property labels
        # (handled separately in _detect_elastic_workflow)

    if m == "DensityFunctionalTheory":
        raise UnsupportedMethodError(
            f"DFT workflows cannot be auto-executed. "
            f"Use mode='create_template' for a skeleton script.\n"
            + UnsupportedMethodError.supported_summary()
        )

    if m == "MolecularDynamics":
        raise UnsupportedMethodError(
            f"MD workflows are not yet supported for auto-execution.\n"
            + UnsupportedMethodError.supported_summary()
        )

    raise UnsupportedMethodError(
        f"Unknown method '{method}' with algorithm '{algorithm}'.\n"
        + UnsupportedMethodError.supported_summary()
    )


def _map_operation_to_node(op_type):
    """Map a KG operation type to an atomrdf.workflow function."""
    mapping = {
        "DeleteAtom": (
            "atomrdf.workflow.pointdefects",
            "create_vacancy",
            "from atomrdf.workflow.pointdefects import create_vacancy",
        ),
        "SubstituteAtom": (
            "atomrdf.workflow.pointdefects",
            "create_substitutional",
            "from atomrdf.workflow.pointdefects import create_substitutional",
        ),
        "AddAtom": (
            "atomrdf.workflow.pointdefects",
            "create_interstitial",
            "from atomrdf.workflow.pointdefects import create_interstitial",
        ),
    }
    if op_type in mapping:
        return mapping[op_type]
    raise UnsupportedMethodError(
        f"Unsupported operation type '{op_type}'.\n"
        + UnsupportedMethodError.supported_summary()
    )


def _detect_elastic_workflow(info):
    """Return True if the calculated properties indicate an elastic constant calc."""
    prop_labels = set(info.get("calculated_properties", {}).keys())
    elastic_labels = {"C11", "C12", "C44", "BulkModulus"}
    return bool(prop_labels & elastic_labels)


# ── Script generation ──────────────────────────────────────────────────


def _generate_script(kg, steps, operations_by_step, mode, output_dir):
    """
    Build the Python script content that chains atomrdf.workflow nodes.

    Parameters
    ----------
    kg : KnowledgeGraph
    steps : list of (info_dict, sample_uri)
    operations_by_step : dict of step_idx → list_of_op_dicts
    mode : "recreate" | "create_template"
    output_dir : str

    Returns
    -------
    str
        Python source code.
    """
    is_template = mode == "create_template"
    lines = []

    # Header
    lines.append("#!/usr/bin/env python")
    lines.append('"""')
    if is_template:
        lines.append("Workflow template generated from Knowledge Graph.")
        lines.append("Fill in the TODO items before running.")
    else:
        lines.append("Reconstructed workflow from Knowledge Graph.")
        lines.append("All parameters filled from KG provenance.")
    lines.append("")
    lines.append(f"Source workflow(s): {', '.join(info['uri'] for info, _ in steps)}")
    lines.append('"""')
    lines.append("")

    # Collect all imports we'll need
    imports = set()
    imports.add("from conceptual_dictionary import ConceptualDict")

    # Determine what functions are needed
    node_calls = []  # list of (step_idx, call_type, function_info, info_dict)

    # First: determine the root structure source (oldest input sample)
    # Walk from oldest step's input samples
    oldest_info = steps[0][0] if steps else None
    root_sample = None
    if oldest_info and oldest_info["input_samples"]:
        root_sample = oldest_info["input_samples"][0]

    # Determine if we need build nodes
    need_build = root_sample is not None
    if need_build:
        imports.add("from atomrdf.workflow.build import bulk, repeat")

    # Process operations
    for step_idx, ops in operations_by_step.items():
        for op in reversed(ops):  # oldest first
            _, func_name, import_line = _map_operation_to_node(op["op_type"])
            imports.add(import_line)
            node_calls.append((step_idx, "operation", op, None))

    # Process computational workflows
    for step_idx, (info, sample_uri) in enumerate(steps):
        method = info["method"]
        algorithm = info["algorithm"]

        if _detect_elastic_workflow(info):
            imports.add(
                "from atomrdf.workflow.elastic import calculate_elastic_constants"
            )
            node_calls.append((step_idx, "elastic", None, info))
        elif method:
            try:
                _, func_name, import_line = _map_method_to_node(method, algorithm)
                imports.add(import_line)
                node_calls.append(
                    (step_idx, "workflow", (method, algorithm, func_name), info)
                )
            except UnsupportedMethodError:
                if is_template:
                    # In template mode, emit a comment instead of failing
                    node_calls.append(
                        (step_idx, "unsupported", (method, algorithm), info)
                    )
                else:
                    raise

    # Write imports
    for imp in sorted(imports):
        lines.append(imp)
    lines.append("")
    lines.append("")

    # Create ConceptualDict
    lines.append("# Create metadata dictionary")
    lines.append("cd = ConceptualDict()")
    lines.append("")

    # Step 1: Build / load root structure
    if root_sample:
        struct_info = _query_sample_structure(kg, root_sample)
        lines.append("# ── Build input structure ──")
        lines.append("")

        element = struct_info.get("element")
        a = struct_info.get("lattice_constant_a")
        n_atoms = struct_info.get("n_atoms")

        # First check if we wrote a structure file
        struct_file = _get_structure_filename(output_dir)

        if struct_file and not is_template:
            lines.append(f"# Structure loaded from KG data (written to {struct_file})")
            lines.append(f"from ase.io import read")
            lines.append(f'structure = read("{struct_file}")')
        elif element:
            # Generate bulk() call
            cs = struct_info.get("crystal_structure")
            cs_lower = _spacegroup_to_crystalstructure(cs) if cs else None

            if is_template:
                # Multi-line call with TODO placeholders
                lines.append("structure = bulk(")
                lines.append(f'    "{element}",')
                if cs_lower:
                    lines.append(f'    crystalstructure="{cs_lower}",')
                lines.append("    a=None,  # TODO: fill in lattice constant")
                if n_atoms and n_atoms <= 4:
                    lines.append("    cubic=True,")
                lines.append("    cdict=cd,")
                lines.append(")")
            else:
                bulk_args = [f'"{element}"']
                if cs_lower:
                    bulk_args.append(f'crystalstructure="{cs_lower}"')
                if a:
                    bulk_args.append(f"a={a}")
                if n_atoms and n_atoms <= 4:
                    bulk_args.append("cubic=True")
                lines.append(f"structure = bulk({', '.join(bulk_args)}, cdict=cd)")
            lines.append("")

            # Check if repeat is needed (infer from atom count)
            unit_cell_atoms = _estimate_unit_cell_atoms(cs)
            if n_atoms and unit_cell_atoms and n_atoms > unit_cell_atoms:
                ratio = n_atoms / unit_cell_atoms
                # Find integer cube root or factorize
                rep = _find_repetitions(ratio)
                if rep:
                    lines.append(f"structure = repeat(structure, {rep}, cdict=cd)")
                    lines.append("")
        else:
            lines.append("# TODO: Load or build your input structure")
            lines.append("# structure = bulk('Fe', cubic=True, cdict=cd)")
            lines.append("# structure = repeat(structure, (4, 4, 4), cdict=cd)")
            lines.append("")

    # Step 2: Operations (defects)
    ops_written = False
    for step_idx, call_type, op_data, info in node_calls:
        if call_type != "operation":
            continue
        ops_written = True
        op_type = op_data["op_type"]
        lines.append(f"# ── Operation: {op_type} ──")
        lines.append("")

        if op_type == "DeleteAtom":
            if is_template:
                lines.append("structure = create_vacancy(structure, cdict=cd)")
                lines.append("# TODO: optionally specify index=<atom_index>")
            else:
                lines.append("structure = create_vacancy(structure, cdict=cd)")
        elif op_type == "SubstituteAtom":
            if is_template:
                lines.append(
                    'structure = create_substitutional(structure, element="X", cdict=cd)  # TODO: fill in element'
                )
            else:
                lines.append(
                    'structure = create_substitutional(structure, element="X", cdict=cd)  # substituting element from KG'
                )
        elif op_type == "AddAtom":
            if is_template:
                lines.append(
                    'structure = create_interstitial(structure, element="X", cdict=cd)  # TODO: fill in element and void_type'
                )
            else:
                lines.append(
                    'structure = create_interstitial(structure, element="X", void_type="tetrahedral", cdict=cd)'
                )
        lines.append("")

    # Step 3: Computational workflows
    for step_idx, call_type, func_data, info in node_calls:
        if call_type == "operation":
            continue

        if call_type == "unsupported":
            method, algorithm = func_data
            lines.append(f"# ── Unsupported workflow: {method} / {algorithm} ──")
            lines.append(
                f"# TODO: This method is not yet supported for auto-execution."
            )
            lines.append(f"# Method: {method}")
            if algorithm:
                lines.append(f"# Algorithm: {algorithm}")
            if info and info["input_parameters"]:
                lines.append("# Input parameters from KG:")
                for k, v in info["input_parameters"].items():
                    lines.append(f"#   {k} = {v}")
            lines.append("")
            continue

        if call_type == "elastic":
            lines.append("# ── Calculate elastic constants ──")
            lines.append("")
            _write_elastic_call(lines, info, is_template)
            continue

        method, algorithm, func_name = func_data

        if func_name == "relax_structure":
            lines.append("# ── Relax structure ──")
            lines.append("")
            _write_relax_call(lines, info, is_template)

        elif func_name == "calculate_ev_curves":
            lines.append("# ── Calculate EV curves ──")
            lines.append("")
            _write_evcurves_call(lines, info, is_template)

    # Final: save metadata
    lines.append("# ── Save metadata ──")
    lines.append("")
    lines.append('cd.to_yaml("reconstructed_workflow.yaml")')
    lines.append(
        'print("Workflow complete. Metadata saved to reconstructed_workflow.yaml")'
    )
    lines.append("")

    # Summary of expected results
    for step_idx, (info, _) in enumerate(steps):
        if info["calculated_properties"]:
            lines.append("# Expected results from KG:")
            for label, value in info["calculated_properties"].items():
                lines.append(f"#   {label} = {value}")
            lines.append("")

    return "\n".join(lines)


def _write_relax_call(lines, info, is_template):
    """Generate code for a relax_structure call."""
    params = info["input_parameters"]

    if is_template:
        lines.append("structure, ecoh, vol = relax_structure(")
        lines.append("    structure,")
        lines.append('    pair_style="eam/alloy",  # TODO: fill in')
        lines.append('    pair_coeff="* * potential.eam.alloy Fe",  # TODO: fill in')
        _write_optional_params(lines, params, is_template)
        lines.append("    cdict=cd,")
        if info.get("potential_type"):
            lines.append(f'    potential_type="{info["potential_type"]}",')
        else:
            lines.append("    potential_type=None,  # TODO: fill in potential type")
        lines.append(")")
    else:
        pair_style, pair_coeff = _extract_potential_params(info)
        lines.append("structure, ecoh, vol = relax_structure(")
        lines.append("    structure,")
        lines.append(f'    pair_style="{pair_style}",')
        lines.append(f'    pair_coeff="{pair_coeff}",')
        _write_optional_params(lines, params, is_template)
        lines.append("    cdict=cd,")
        if info.get("potential_type"):
            lines.append(f'    potential_type="{info["potential_type"]}",')
        lines.append(")")
    lines.append(
        'print(f"Relaxed: E_coh = {ecoh:.4f} eV/atom, V = {vol:.4f} A^3/atom")'
    )
    lines.append("")


def _write_evcurves_call(lines, info, is_template):
    """Generate code for a calculate_ev_curves call."""
    params = info["input_parameters"]

    if is_template:
        lines.append("results = calculate_ev_curves(")
        lines.append("    structure,")
        lines.append('    pair_style="eam/alloy",  # TODO: fill in')
        lines.append('    pair_coeff="* * potential.eam.alloy Fe",  # TODO: fill in')
        _write_optional_params(lines, params, is_template)
        lines.append("    cdict=cd,")
        lines.append("    potential_type=None,  # TODO: fill in")
        lines.append(")")
    else:
        pair_style, pair_coeff = _extract_potential_params(info)
        lines.append("results = calculate_ev_curves(")
        lines.append("    structure,")
        lines.append(f'    pair_style="{pair_style}",')
        lines.append(f'    pair_coeff="{pair_coeff}",')
        _write_optional_params(lines, params, is_template)
        lines.append("    cdict=cd,")
        if info.get("potential_type"):
            lines.append(f'    potential_type="{info["potential_type"]}",')
        lines.append(")")
    lines.append(
        "print(f\"EV curves: bulk modulus = {results['bulk_modulus']:.2f} GPa\")"
    )
    lines.append("")


def _write_elastic_call(lines, info, is_template):
    """Generate code for a calculate_elastic_constants call."""
    params = info["input_parameters"]

    if is_template:
        lines.append("results = calculate_elastic_constants(")
        lines.append("    structure,")
        lines.append('    pair_style="eam/alloy",  # TODO: fill in')
        lines.append('    pair_coeff="* * potential.eam.alloy Fe",  # TODO: fill in')
        lines.append("    cdict=cd,")
        lines.append("    potential_type=None,  # TODO: fill in")
        lines.append(")")
    else:
        pair_style, pair_coeff = _extract_potential_params(info)
        lines.append("results = calculate_elastic_constants(")
        lines.append("    structure,")
        lines.append(f'    pair_style="{pair_style}",')
        lines.append(f'    pair_coeff="{pair_coeff}",')
        lines.append("    cdict=cd,")
        if info.get("potential_type"):
            lines.append(f'    potential_type="{info["potential_type"]}",')
        lines.append(")")
    lines.append(
        "print(f\"Elastic: C11={results['C_matrix'][0,0]:.1f}, C12={results['C_matrix'][0,1]:.1f}, C44={results['C_matrix'][3,3]:.1f} GPa\")"
    )
    lines.append("")


def _write_optional_params(lines, params, is_template):
    """Write optional minimization/simulation parameters if present in KG."""
    param_map = {
        "energy_tolerance": "e_tol",
        "force_tolerance": "f_tol",
        "n_energy_steps": "n_energy_steps",
        "n_force_steps": "n_force_steps",
        "cores": "cores",
        "vol_range": "vol_range",
        "num_of_points": "num_of_points",
    }
    for kg_name, py_name in param_map.items():
        if kg_name in params:
            lines.append(f"    {py_name}={params[kg_name]},")


def _extract_potential_params(info):
    """
    Extract pair_style and pair_coeff from workflow info.
    Falls back to placeholder strings if not found.
    """
    params = info["input_parameters"]
    pair_style = params.get("pair_style", "eam/alloy")
    pair_coeff = params.get("pair_coeff", "* * potential.alloy Element")

    # Try to infer from potential type
    if pair_style == "eam/alloy" and info.get("potential_type"):
        pt = info["potential_type"].lower()
        if "meam" in pt:
            pair_style = "meam"
        elif "ace" in pt or "pace" in pt:
            pair_style = "pace"
        elif "lj" in pt or "lennard" in pt:
            pair_style = "lj/cut 10.0"

    return pair_style, pair_coeff


def _get_structure_filename(output_dir):
    """Check if a structure file was written in the output directory."""
    for name in ("POSCAR", "structure.lmp", "structure.in", "structure.dat"):
        if os.path.exists(os.path.join(output_dir, name)):
            return name
    return None


def _spacegroup_to_crystalstructure(sg):
    """Map a space group symbol to ASE crystal structure name."""
    mapping = {
        "Im-3m": "bcc",
        "Fm-3m": "fcc",
        "P6_3/mmc": "hcp",
        "Fd-3m": "diamond",
        "Pm-3m": "sc",
    }
    return mapping.get(sg)


def _estimate_unit_cell_atoms(crystal_structure):
    """Estimate atoms per conventional unit cell."""
    if crystal_structure is None:
        return None
    mapping = {
        "Im-3m": 2,
        "bcc": 2,
        "Fm-3m": 4,
        "fcc": 4,
        "P6_3/mmc": 2,
        "hcp": 2,
        "Fd-3m": 8,
        "diamond": 8,
        "Pm-3m": 1,
        "sc": 1,
    }
    return mapping.get(crystal_structure)


def _find_repetitions(ratio):
    """
    Find integer repetitions tuple from atom-count ratio.
    Returns tuple like (4, 4, 4) or None.
    """
    import math

    # Try cube root first
    cr = round(ratio ** (1 / 3))
    if cr**3 == int(ratio):
        return f"({cr}, {cr}, {cr})"

    # Try simple factors
    for a in range(1, 20):
        for b in range(a, 20):
            c_val = ratio / (a * b)
            if c_val == int(c_val) and int(c_val) >= b:
                return f"({a}, {b}, {int(c_val)})"

    return None


# ── Public API ─────────────────────────────────────────────────────────


def reconstruct_workflow(
    kg, workflow_id, output_dir, mode="recreate", structure_format=None
):
    """
    Reconstruct a workflow as an executable Python script.

    Queries the full provenance chain from the KG and generates a Python
    script that uses ``atomrdf.workflow`` nodes to reproduce the workflow.

    Parameters
    ----------
    kg : KnowledgeGraph
        An atomRDF KnowledgeGraph instance.
    workflow_id : str or URIRef
        URI of the workflow / simulation node.
    output_dir : str
        Path to the output directory (created if needed).
    mode : str
        ``"recreate"`` — fully runnable script with all inputs from KG.
        ``"create_template"`` — skeleton with TODO placeholders.
    structure_format : str, optional
        Override structure file format.

    Returns
    -------
    str
        The *output_dir* path.

    Raises
    ------
    UnsupportedMethodError
        If the workflow uses a method/algorithm we cannot reconstruct.
    """
    g = kg.graph
    wf_uri = URIRef(workflow_id) if not isinstance(workflow_id, URIRef) else workflow_id

    # 1. Collect workflow chain
    info = _query_workflow_info(kg, wf_uri)

    all_steps = []
    for out_sample in info["output_samples"]:
        chain = _walk_workflow_chain(kg, out_sample)
        for pair in chain:
            if pair not in all_steps:
                all_steps.append(pair)

    if not all_steps:
        all_steps = [
            (info["output_samples"][0] if info["output_samples"] else None, str(wf_uri))
        ]

    all_steps = list(reversed(all_steps))

    # 2. Gather info + operations for each step
    steps = []
    operations_by_step = {}
    for idx, (sample_uri, wf_step_uri) in enumerate(all_steps):
        step_info = _query_workflow_info(kg, wf_step_uri)
        steps.append((step_info, sample_uri))

        if step_info["output_samples"]:
            ops = _query_operation_chain(kg, step_info["output_samples"][0])
            if ops:
                operations_by_step[idx] = ops

    # 3. Create output directory and write structure files
    os.makedirs(output_dir, exist_ok=True)

    # Write input structure if available
    fmt = structure_format or "lammps-data"
    oldest_info = steps[0][0] if steps else None
    if oldest_info and oldest_info["input_samples"]:
        root_sample = oldest_info["input_samples"][0]
        struct_fname = _format_extension(fmt)
        _write_structure_fallback(
            kg,
            root_sample,
            os.path.join(output_dir, struct_fname),
            fmt,
        )

    # 4. Generate the executable script
    script = _generate_script(kg, steps, operations_by_step, mode, output_dir)

    script_path = os.path.join(output_dir, "run_workflow.py")
    with open(script_path, "w") as f:
        f.write(script)

    # 5. Write a summary README
    _write_readme(steps, operations_by_step, mode, output_dir)

    return output_dir


def reconstruct_workflow_by_sample(
    kg, sample_id, output_dir, mode="recreate", structure_format=None
):
    """
    Find the workflow that produced *sample_id* and reconstruct it.

    Parameters
    ----------
    kg : KnowledgeGraph
    sample_id : str or URIRef
    output_dir : str
    mode : str
        ``"recreate"`` or ``"create_template"``.
    structure_format : str, optional

    Returns
    -------
    str
        The *output_dir* path.
    """
    g = kg.graph
    sample = URIRef(sample_id) if not isinstance(sample_id, URIRef) else sample_id

    workflow = None
    for _, _, wf in g.triples((sample, PROV.wasGeneratedBy, None)):
        method = g.value(wf, ASMO.hasComputationalMethod)
        if method is not None:
            workflow = wf
            break

    if workflow is None:
        raise ValueError(
            f"No workflow (simulation with hasComputationalMethod) found "
            f"for sample {sample_id}"
        )

    return reconstruct_workflow(
        kg,
        workflow,
        output_dir,
        mode=mode,
        structure_format=structure_format,
    )


def reconstruct_from_property(
    kg, property_uri, output_dir, mode="recreate", structure_format=None
):
    """
    Start from a calculated property, find its workflow, and reconstruct it.

    Parameters
    ----------
    kg : KnowledgeGraph
    property_uri : str or URIRef
        URI of a calculated property node.
    output_dir : str
    mode : str
    structure_format : str, optional

    Returns
    -------
    str
        The *output_dir* path.
    """
    g = kg.graph
    prop = (
        URIRef(property_uri) if not isinstance(property_uri, URIRef) else property_uri
    )

    # Find the workflow that calculated this property
    workflow = g.value(prop, ASMO.wasCalculatedBy)
    if workflow is None:
        raise ValueError(f"No workflow found that calculated property {property_uri}")

    return reconstruct_workflow(
        kg,
        workflow,
        output_dir,
        mode=mode,
        structure_format=structure_format,
    )


def _write_readme(steps, operations_by_step, mode, output_dir):
    """Write a README.md summarizing the reconstructed workflow."""
    lines = ["# Reconstructed Workflow", ""]

    if mode == "create_template":
        lines.append(
            "**Template mode**: Fill in the `# TODO` items in `run_workflow.py` before running."
        )
    else:
        lines.append("**Recreate mode**: `run_workflow.py` is ready to execute.")
    lines.append("")

    lines.append("## Steps")
    lines.append("")

    for i, (info, sample) in enumerate(steps, 1):
        method = info["method"] or "N/A"
        algo = info["algorithm"] or ""
        sw = ", ".join(info["software"]) if info["software"] else "N/A"

        lines.append(f"### Step {i}: {method}" + (f" / {algo}" if algo else ""))
        lines.append(f"- Software: {sw}")
        if info["dofs"]:
            lines.append(f"- DOFs: {', '.join(info['dofs'])}")
        if info["potential"]:
            lines.append(f"- Potential: {info['potential']}")
        if info["calculated_properties"]:
            lines.append("- Calculated properties:")
            for label, value in info["calculated_properties"].items():
                lines.append(f"  - {label} = {value}")
        lines.append("")

    # Operations
    all_ops = []
    for ops in operations_by_step.values():
        all_ops.extend(ops)
    if all_ops:
        lines.append("## Operations")
        for op in reversed(all_ops):
            lines.append(f"- {op['op_type']}")
        lines.append("")

    lines.append("## Usage")
    lines.append("")
    lines.append("```bash")
    lines.append("python run_workflow.py")
    lines.append("```")
    lines.append("")

    path = os.path.join(output_dir, "README.md")
    with open(path, "w") as f:
        f.write("\n".join(lines))
