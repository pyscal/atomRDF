"""Code generation from atomRDF provenance chains.

Translates a :class:`~atomrdf.io.provenance.Provenance` object into an
executable Python script that reconstructs the workflow using real
workflow functions (e.g. ``atomrdf.workflow.lammps``).
"""

import os
import re
from collections import OrderedDict

from rdflib import URIRef, RDFS

from atomrdf.namespace import ASMO


# ------------------------------------------------------------------ #
# Helpers                                                              #
# ------------------------------------------------------------------ #


def _sanitize(name):
    """Convert *name* to a valid, snake_case Python identifier."""
    # Insert underscore at camelCase boundaries
    name = re.sub(r"(?<=[a-z])([A-Z])", r"_\1", str(name))
    name = re.sub(r"[^A-Za-z0-9_]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    if not name or name[0].isdigit():
        name = "v_" + name
    return name.lower()


def _sample_label(kg, sample_id):
    """Human-readable label for a sample URI."""
    label = kg.get_label(URIRef(str(sample_id)))
    if label:
        return str(label)
    s = str(sample_id)
    if ":" in s:
        return s.split(":")[-1][:8]
    return s.rsplit("/", 1)[-1][:8]


# ------------------------------------------------------------------ #
# CodeContext                                                          #
# ------------------------------------------------------------------ #


class CodeContext:
    """Accumulates generated Python code and tracks variable mappings."""

    def __init__(self, kg):
        self.kg = kg
        self.uri_to_var: dict[str, str] = {}
        self.user_inputs: OrderedDict[str, str] = OrderedDict()
        self._imports: list[str] = []
        self.lines: list[str] = []
        self.structures: dict[str, object] = {}  # filename → ASE Atoms
        self._var_counts: dict[str, int] = {}
        self._step_num = 0

    # -- variable naming --------------------------------------------- #

    def make_var(self, hint, uri=None):
        """Create a unique Python variable name from *hint*."""
        name = _sanitize(hint)
        count = self._var_counts.get(name, 0) + 1
        self._var_counts[name] = count
        full = f"{name}_{count}" if count > 1 else name
        if uri:
            self.uri_to_var[str(uri)] = full
        return full

    # -- accumulation helpers ---------------------------------------- #

    def add_user_input(self, name, description):
        if name not in self.user_inputs:
            self.user_inputs[name] = description

    def add_import(self, line):
        if line not in self._imports:
            self._imports.append(line)

    def comment(self, text):
        self._step_num += 1
        self.lines.append(f"\n# === Step {self._step_num}: {text} ===")

    def code(self, line):
        self.lines.append(line)

    def add_structure(self, filename, atoms):
        self.structures[filename] = atoms

    # -- rendering --------------------------------------------------- #

    def render(self):
        """Assemble the full workflow script as a string."""
        parts = [
            '"""',
            "Auto-generated workflow script from atomRDF provenance.",
            '"""',
            "",
        ]

        if self.user_inputs:
            parts.append("# === User-fillable parameters ===")
            for name, desc in self.user_inputs.items():
                parts.append(f'{name} = "..."  # {desc}')
            parts.append("")

        if self._imports:
            for imp in self._imports:
                parts.append(imp)
            parts.append("")

        parts.extend(self.lines)
        parts.append("")

        return "\n".join(parts)

    def write(self, output_dir, script=None):
        """Write workflow script and structure files to *output_dir*."""
        os.makedirs(output_dir, exist_ok=True)

        if script is None:
            script = self.render()
        with open(os.path.join(output_dir, "workflow.py"), "w") as f:
            f.write(script)

        if self.structures:
            struct_dir = os.path.join(output_dir, "structures")
            os.makedirs(struct_dir, exist_ok=True)
            from ase.io import write as ase_write

            for fname, atoms in self.structures.items():
                ase_write(os.path.join(struct_dir, fname), atoms, format="json")


# ------------------------------------------------------------------ #
# Code generation entry point                                         #
# ------------------------------------------------------------------ #


def generate_code(provenance, output_dir=None):
    """Generate an executable Python script from a Provenance object.

    Parameters
    ----------
    provenance : Provenance
    output_dir : str, optional
        If given, write ``workflow.py`` and ``structures/`` here.

    Returns
    -------
    str
        The generated Python script.
    """
    ctx = CodeContext(provenance.kg)

    # 1. Process chain steps (sample → activity → sample)
    chain_steps = [s for s in provenance._steps if "result_property" not in s]
    for step in chain_steps:
        _handle_chain_step(provenance, ctx, step)

    # 2. Ensure sub-chains needed by math operands are generated
    math_steps = [s for s in provenance._steps if "result_property" in s]
    _ensure_sub_chains(provenance, ctx, math_steps)

    # 3. Generate arithmetic for math operations
    for step in math_steps:
        _handle_math_step(provenance, ctx, step)

    script = ctx.render()
    if output_dir is not None:
        ctx.write(output_dir, script)
    return script


# ------------------------------------------------------------------ #
# Step handlers                                                        #
# ------------------------------------------------------------------ #


def _handle_chain_step(provenance, ctx, step):
    """Dispatch a single chain step to the appropriate handler."""
    out_id = step["output_sample_id"]
    if out_id in ctx.uri_to_var:
        return  # already processed

    in_id = step.get("input_sample_id")
    if in_id and in_id not in ctx.uri_to_var:
        _ensure_structure(ctx, in_id, step.get("input_sample"))

    act_type = step["activity_type"]
    if act_type == "Simulation":
        _handle_simulation(provenance, ctx, step)
    elif act_type in (
        "SubstituteAtom",
        "DeleteAtom",
        "AddAtom",
        "Rotate",
        "Translate",
        "Shear",
    ):
        _handle_operation(ctx, step)
    else:
        _ensure_structure(ctx, out_id, step.get("output_sample"))


def _ensure_structure(ctx, sample_id, atoms):
    """Write a structure file and emit a ``read()`` call."""
    if str(sample_id) in ctx.uri_to_var:
        return

    label = _sample_label(ctx.kg, sample_id)
    var = ctx.make_var(f"atoms_{label}", uri=sample_id)
    fname = f"{_sanitize(label)}.json"

    ctx.add_import("from ase.io import read")
    ctx.comment(f"Load structure: {label}")
    ctx.code(f'{var} = read("structures/{fname}")')

    if atoms is not None:
        ctx.add_structure(fname, atoms)
    else:
        ctx.code(
            f"# WARNING: structure not available; provide structures/{fname} manually"
        )


def _handle_simulation(provenance, ctx, step):
    """Generate code for a Simulation step."""
    in_id = step.get("input_sample_id")
    out_id = step["output_sample_id"]
    in_var = ctx.uri_to_var.get(in_id) if in_id else None

    if in_var is None:
        # Root step — no input, just load the output structure
        _ensure_structure(ctx, out_id, step.get("output_sample"))
        _register_calc_properties(provenance, ctx, step, var_map=None)
        return

    software = step.get("software", [])
    method = step.get("method") or "Simulation"
    dof = step.get("degrees_of_freedom", [])
    out_label = _sample_label(ctx.kg, out_id)

    is_lammps = any("LAMMPS" in str(s) for s in software)
    has_cell_relax = any("CellVolume" in str(d) for d in dof)

    if is_lammps:
        ctx.add_user_input("pair_style", "LAMMPS pair style")
        ctx.add_user_input("pair_coeff", "LAMMPS pair coefficients")

        if has_cell_relax:
            func = "calculate_energy_relax"
            ctx.add_import("from atomrdf.workflow.lammps import calculate_energy_relax")
            out_var = ctx.make_var(f"atoms_{out_label}", uri=out_id)
            ecoh_var = ctx.make_var("ecoh")
            vol_var = ctx.make_var("vol")

            ctx.comment(f"{method} \u2192 {out_label} (relax)")
            ctx.code(f"{out_var}, {ecoh_var}, {vol_var} = {func}(")
            ctx.code(f"    {in_var}, pair_style=pair_style, pair_coeff=pair_coeff")
            ctx.code(")")

        else:
            func = "calculate_energy_rigid"
            ctx.add_import("from atomrdf.workflow.lammps import calculate_energy_rigid")
            ecoh_var = ctx.make_var("ecoh")
            vol_var = ctx.make_var("vol")

            ctx.comment(f"{method} \u2192 {out_label} (rigid)")
            ctx.code(f"{ecoh_var}, {vol_var} = {func}(")
            ctx.code(f"    {in_var}, pair_style=pair_style, pair_coeff=pair_coeff")
            ctx.code(")")

            # Rigid doesn't return a new structure; reuse input variable
            ctx.uri_to_var[out_id] = in_var

        _register_calc_properties(
            provenance,
            ctx,
            step,
            var_map={"energy": ecoh_var, "volume": vol_var},
        )

    else:
        # Unknown simulation software — emit placeholder
        sw_names = ", ".join(software) if software else "unknown"
        ctx.comment(f"{method} ({sw_names})")
        ctx.code(f"# TODO: implement simulation for {sw_names}")
        _ensure_structure(ctx, out_id, step.get("output_sample"))
        _register_calc_properties(provenance, ctx, step, var_map=None)


def _register_calc_properties(provenance, ctx, step, var_map):
    """Map calculated-property URIs to the variables holding their values."""
    out_id = step["output_sample_id"]
    sim_id = step["activity_id"]

    for _, _, prop_uri in provenance.kg.triples(
        (URIRef(out_id), ASMO.hasCalculatedProperty, None)
    ):
        calc_by = provenance.kg.graph.value(prop_uri, ASMO.wasCalculatedBy)
        if calc_by is None or str(calc_by) != sim_id:
            continue

        label_node = provenance.kg.graph.value(prop_uri, RDFS.label)
        label_str = str(label_node) if label_node else ""

        if var_map is not None:
            if any(k in label_str.lower() for k in ("energy",)):
                ctx.uri_to_var[str(prop_uri)] = var_map["energy"]
            elif any(k in label_str.lower() for k in ("volume",)):
                ctx.uri_to_var[str(prop_uri)] = var_map["volume"]
        else:
            # No simulation call — use stored value as constant
            val_node = provenance.kg.graph.value(prop_uri, ASMO.hasValue)
            if val_node is not None:
                var = ctx.make_var(_sanitize(label_str or "prop"), uri=prop_uri)
                ctx.code(f"{var} = {float(val_node)}  # stored from graph")


def _handle_operation(ctx, step):
    """Generate code for an atom operation (substitute, delete, etc.)."""
    out_id = step["output_sample_id"]
    act_type = step["activity_type"]
    out_label = _sample_label(ctx.kg, out_id)

    var = ctx.make_var(f"atoms_{out_label}", uri=out_id)
    fname = f"{_sanitize(out_label)}.json"

    ctx.add_import("from ase.io import read")
    ctx.comment(f"{act_type} \u2192 {out_label}")
    ctx.code(f'{var} = read("structures/{fname}")')

    atoms = step.get("output_sample")
    if atoms is not None:
        ctx.add_structure(fname, atoms)
    else:
        ctx.code(
            f"# WARNING: structure not available; provide structures/{fname} manually"
        )


# ------------------------------------------------------------------ #
# Math operations                                                      #
# ------------------------------------------------------------------ #


def _handle_math_step(provenance, ctx, step):
    """Generate arithmetic code for a math operation."""
    activity = step.get("activity")
    result = step.get("result_property", {})
    act_type = step["activity_type"]

    if activity is None:
        return

    result_label = result.get("label", "result")
    result_var = ctx.make_var(_sanitize(result_label), uri=result.get("uri"))

    def resolve(operand):
        if operand is None:
            return "None  # missing operand"
        if isinstance(operand, (int, float)):
            return repr(operand)
        key = str(operand)
        if key in ctx.uri_to_var:
            return ctx.uri_to_var[key]
        # Fallback: use stored value from graph
        val_node = provenance.kg.graph.value(URIRef(key), ASMO.hasValue)
        if val_node is not None:
            return repr(float(val_node))
        return f"UNKNOWN  # could not resolve: {key}"

    ctx.comment(act_type)

    if act_type == "Subtraction":
        a = resolve(getattr(activity, "minuend", None))
        b = resolve(getattr(activity, "subtrahend", None))
        ctx.code(f"{result_var} = {a} - {b}")

    elif act_type == "Addition":
        addends = getattr(activity, "addend", []) or []
        parts = [resolve(a) for a in addends]
        ctx.code(f"{result_var} = {' + '.join(parts) if parts else '0'}")

    elif act_type == "Multiplication":
        factors = getattr(activity, "factor", []) or []
        parts = [resolve(f) for f in factors]
        ctx.code(f"{result_var} = {' * '.join(parts) if parts else '1'}")

    elif act_type == "Division":
        a = resolve(getattr(activity, "dividend", None))
        b = resolve(getattr(activity, "divisor", None))
        ctx.code(f"{result_var} = {a} / {b}")

    elif act_type == "Exponentiation":
        base_val = resolve(getattr(activity, "base", None))
        exp_val = resolve(getattr(activity, "exponent", None))
        ctx.code(f"{result_var} = {base_val} ** {exp_val}")

    else:
        ctx.code(f"# TODO: unsupported math operation: {act_type}")

    # Annotate with expected result
    if result.get("value") is not None:
        unit = result.get("unit", "") or ""
        ctx.code(f"# Expected: {result['value']:.6g} {unit}".rstrip())


# ------------------------------------------------------------------ #
# Sub-chain handling                                                   #
# ------------------------------------------------------------------ #


def _ensure_sub_chains(provenance, ctx, math_steps):
    """Trace sub-chains for operand properties owned by external samples."""
    from atomrdf.io.provenance import Provenance

    known_ids = set()
    for s in provenance._steps:
        if s.get("output_sample_id"):
            known_ids.add(s["output_sample_id"])
        if s.get("input_sample_id"):
            known_ids.add(s["input_sample_id"])

    for step in math_steps:
        activity = step.get("activity")
        if activity is None:
            continue
        for operand in _collect_operands(activity):
            if not isinstance(operand, str) or operand in ctx.uri_to_var:
                continue
            owner = provenance.kg.graph.value(
                None, ASMO.hasCalculatedProperty, URIRef(operand)
            )
            if owner is None or str(owner) in known_ids:
                continue
            try:
                sub_prov = Provenance.from_sample(provenance.kg, owner)
                for sub_step in sub_prov._steps:
                    if "result_property" not in sub_step:
                        _handle_chain_step(provenance, ctx, sub_step)
            except Exception:
                pass


def _collect_operands(activity):
    """Gather all operand values from a math-operation activity object."""
    result = []
    for attr in (
        "minuend",
        "subtrahend",
        "addend",
        "factor",
        "dividend",
        "divisor",
        "base",
        "exponent",
    ):
        val = getattr(activity, attr, None)
        if val is None:
            continue
        if isinstance(val, list):
            result.extend(val)
        else:
            result.append(val)
    return result
