"""Code generation from atomRDF provenance chains.

Translates a :class:`~atomrdf.io.provenance.Provenance` object into an
executable Python script that reconstructs the workflow using real
workflow functions (e.g. ``atomrdf.workflow.lammps``).

Simulation dispatch mapping
---------------------------
Code selection in ``_handle_simulation()`` is driven by three ASMO-mapped
attributes read from the RDF graph.  The table below shows every recognised
combination and the generated output.  Rows marked ``# TODO`` have no
implementing function yet.

Software (step["software"])
    Detected by: ``any("<NAME>" in s for s in software)``
    ASMO class  : ``PROV.SoftwareAgent``  (label stored as RDFS.label)

Method (step["method"])
    ASMO classes: ``ASMO.MolecularStatics``
                  ``ASMO.MolecularDynamics``
                  ``ASMO.DensityFunctionalTheory``

Degrees of freedom (step["degrees_of_freedom"])
    ASMO named individuals:
        ``ASMO.AtomicPositionRelaxation``   basename="AtomicPositionRelaxation"
        ``ASMO.CellVolumeRelaxation``       basename="CellVolumeRelaxation"
        ``ASMO.CellShapeRelaxation``        basename="CellShapeRelaxation"

Algorithm (step["algorithm"])
    ASMO classes: ``ASMO.EquationOfStateFit``
                  ``ASMO.QuasiHarmonicApproximation``
                  ``ASMO.ThermodynamicIntegration``
                  ``ASMO.ANNNImodel``
                  ``ASMO.TensileTest``
                  ``ASMO.CompressionTest``

Thermodynamic ensemble (step["activity"].thermodynamic_ensemble)
    ASMO named individuals:
        ``ASMO.CanonicalEnsemble``              (NVT)
        ``ASMO.MicrocanonicalEnsemble``         (NVE)
        ``ASMO.IsothermalIsobaricEnsemble``     (NPT)
        ``ASMO.IsoenthalpicIsobaricEnsemble``   (NPH)
        ``ASMO.GrandCanonicalEnsemble``         (muVT)

Interatomic potential (step["interatomic_potential"]["type"])
    ASMO classes: ``ASMO.InteratomicPotential``      (generic)
                  ``ASMO.EmbeddedAtomModel``          (EAM)
                  ``ASMO.ModifiedEmbeddedAtomModel``  (MEAM)
                  ``ASMO.LennardJonesPotential``      (LJ)
                  ``ASMO.MachineLearningPotential``   (MLP / GRACE / ACE etc.)

XC functional  (step["activity"].xc_functional) — DFT only
    MDO classes:  ``MDO.ExchangeCorrelationEnergyFunctional``  (generic)
                  ``MDO.GeneralizedGradientApproximation``     (GGA)
                  ``MDO.LocalDensityApproximation``            (LDA)
                  ``MDO.HybridFunctional``
                  ``MDO.HybridGeneralizedGradientApproximation``
                  ``MDO.HybridMetaGeneralizedGradientApproximation``
                  ``MDO.MetaGeneralizedGradientApproximation``

Dispatch table
--------------
Software   Method               DOFs                    Algorithm              Generated call
---------- -------------------- ----------------------- ---------------------- -------------------------------------------------------
LAMMPS     MolecularStatics     AtomicPosition          —                      calculate_energy_rigid(atoms, ...)             ✓
LAMMPS     MolecularStatics     AtomicPosition+CellVol  —                      calculate_energy_relax(atoms, ...)             ✓
LAMMPS     MolecularDynamics    AtomicPosition          NVT/NVE/NPT/NPH        # TODO: calculate_md(atoms, ensemble, ...)
LAMMPS     MolecularStatics     —                       EquationOfStateFit     # TODO: calculate_eos(atoms, ...)
LAMMPS     MolecularDynamics    AtomicPosition          ThermodynamicInteg.    # TODO: calculate_ti(atoms, ...)
LAMMPS     MolecularDynamics    AtomicPosition          QuasiHarmonicApprox.   # TODO: calculate_qha(atoms, ...)
LAMMPS     MolecularStatics     —                       TensileTest            # TODO: calculate_tensile(atoms, ...)
LAMMPS     MolecularStatics     —                       CompressionTest        # TODO: calculate_compression(atoms, ...)
VASP       DFT                  AtomicPosition          —                      # TODO: run_vasp(atoms, xc, ...)
VASP       DFT                  AtomicPosition+CellVol  —                      # TODO: run_vasp(atoms, xc, relax=True, ...)
QuantumESP DFT                  AtomicPosition          —                      # TODO: run_qe(atoms, xc, ...)
(other)    (any)                (any)                   (any)                  # TODO placeholder emitted
"""

import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Optional

from rdflib import URIRef, RDFS

from atomrdf.namespace import ASMO


# ------------------------------------------------------------------ #
# Simulation dispatch table                                            #
# ------------------------------------------------------------------ #


@dataclass
class WorkflowNode:
    """Describes how to generate code for one class of simulation step.

    Matching rules (all non-None fields must match for the handler to fire):
      software   – substring found in any element of step["software"]
      method     – exact basename of step["method"]
      algorithm  – exact basename of step["algorithm"]
      dof        – frozenset of DOF basenames that must ALL be present
                   in step["degrees_of_freedom"]

    Generation fields:
      import_line       – import statement to add (None → emit # TODO)
      func              – callable name (None → emit # TODO)
      returns_structure – True  → ``out, ecoh, vol = func(in, ...)``
                          False → ``ecoh, vol = func(in, ...)``
      user_inputs       – {param_name: description} added to script header
      call_kwargs       – extra keyword-args string passed after the atoms arg
    """

    note: str
    software: Optional[str] = None  # substring to match in software labels
    method: Optional[str] = None  # method basename
    algorithm: Optional[str] = None  # algorithm basename
    dof: Optional[frozenset] = None  # DOF basenames required (all must match)
    import_line: Optional[str] = None  # None → TODO placeholder
    func: Optional[str] = None  # None → TODO placeholder
    returns_structure: bool = True
    user_inputs: dict = field(default_factory=dict)
    call_kwargs: str = ""


# Entries are checked in order; the first match wins.
# More-specific entries (with DOF / algorithm constraints) must come first.
_DISPATCH_TABLE = [
    # ---------------------------------------------------------------- #
    # LAMMPS  — algorithm-specific entries MUST come before generic    #
    # ones, because _match_handler() returns the first match.          #
    # ---------------------------------------------------------------- #
    # --- EquationOfStateFit -----------------------------------------
    WorkflowNode(
        note="LAMMPS equation-of-state fit (Birch-Murnaghan)",
        software="LAMMPS",
        method="MolecularStatics",
        algorithm="EquationOfStateFit",
        import_line="from workflows.evcurves import calculate_ev_curves",
        func="calculate_ev_curves",
        returns_structure=True,
        user_inputs={
            "pair_style": "LAMMPS pair style",
            "pair_coeff": "LAMMPS pair coefficients",
        },
        call_kwargs="pair_style=pair_style, pair_coeff=pair_coeff",
    ),
    # --- QuasiHarmonicApproximation ---------------------------------
    WorkflowNode(
        note="LAMMPS quasi-harmonic approximation",
        software="LAMMPS",
        algorithm="QuasiHarmonicApproximation",
        # TODO: func="calculate_qha", import_line="from atomrdf.workflow.lammps import calculate_qha"
    ),
    # --- ThermodynamicIntegration ------------------------------------
    WorkflowNode(
        note="LAMMPS thermodynamic integration",
        software="LAMMPS",
        algorithm="ThermodynamicIntegration",
        # TODO: func="calculate_ti", import_line="from atomrdf.workflow.lammps import calculate_ti"
    ),
    # --- TensileTest --------------------------------------------------
    WorkflowNode(
        note="LAMMPS tensile test",
        software="LAMMPS",
        algorithm="TensileTest",
        # TODO: func="calculate_tensile", import_line="from atomrdf.workflow.lammps import calculate_tensile"
    ),
    # --- CompressionTest ---------------------------------------------
    WorkflowNode(
        note="LAMMPS compression test",
        software="LAMMPS",
        algorithm="CompressionTest",
        # TODO: func="calculate_compression", import_line="from atomrdf.workflow.lammps import calculate_compression"
    ),
    # --- MolecularDynamics (generic, catch-all ensemble) ------------
    WorkflowNode(
        note="LAMMPS MolecularDynamics (NVT/NVE/NPT/NPH — pick ensemble from activity)",
        software="LAMMPS",
        method="MolecularDynamics",
        # TODO: func="calculate_md", import_line="from atomrdf.workflow.lammps import calculate_md"
    ),
    # --- MolecularStatics with cell relaxation ----------------------
    WorkflowNode(
        note="LAMMPS MolecularStatics with cell relaxation",
        software="LAMMPS",
        method="MolecularStatics",
        dof=frozenset({"CellVolumeRelaxation"}),
        import_line="from atomrdf.workflow.lammps import calculate_energy_relax",
        func="calculate_energy_relax",
        returns_structure=True,
        user_inputs={
            "pair_style": "LAMMPS pair style",
            "pair_coeff": "LAMMPS pair coefficients",
        },
        call_kwargs="pair_style=pair_style, pair_coeff=pair_coeff",
    ),
    # --- MolecularStatics rigid (catch-all, no DOF/algorithm) --------
    WorkflowNode(
        note="LAMMPS MolecularStatics rigid (atoms-only relaxation)",
        software="LAMMPS",
        method="MolecularStatics",
        import_line="from atomrdf.workflow.lammps import calculate_energy_rigid",
        func="calculate_energy_rigid",
        returns_structure=False,
        user_inputs={
            "pair_style": "LAMMPS pair style",
            "pair_coeff": "LAMMPS pair coefficients",
        },
        call_kwargs="pair_style=pair_style, pair_coeff=pair_coeff",
    ),
    # ---------------------------------------------------------------- #
    # VASP                                                              #
    # ---------------------------------------------------------------- #
    WorkflowNode(
        note="VASP DFT with cell relaxation",
        software="VASP",
        method="DensityFunctionalTheory",
        dof=frozenset({"CellVolumeRelaxation"}),
        # TODO: func="run_vasp", import_line="from atomrdf.workflow.vasp import run_vasp"
    ),
    WorkflowNode(
        note="VASP DFT single-point / ionic relaxation",
        software="VASP",
        method="DensityFunctionalTheory",
        # TODO: func="run_vasp", import_line="from atomrdf.workflow.vasp import run_vasp"
    ),
    # ---------------------------------------------------------------- #
    # Quantum ESPRESSO                                                  #
    # ---------------------------------------------------------------- #
    WorkflowNode(
        note="Quantum ESPRESSO DFT with cell relaxation",
        software="QuantumESPRESSO",
        method="DensityFunctionalTheory",
        dof=frozenset({"CellVolumeRelaxation"}),
        # TODO: func="run_qe", import_line="from atomrdf.workflow.qe import run_qe"
    ),
    WorkflowNode(
        note="Quantum ESPRESSO DFT single-point / ionic relaxation",
        software="QuantumESPRESSO",
        method="DensityFunctionalTheory",
        # TODO: func="run_qe", import_line="from atomrdf.workflow.qe import run_qe"
    ),
]


def _match_handler(step) -> Optional["WorkflowNode"]:
    """Return the first entry in _DISPATCH_TABLE that matches *step*."""
    software = step.get("software", [])
    method = step.get("method")
    algorithm = step.get("algorithm")
    dof = set(step.get("degrees_of_freedom", []))

    for handler in _DISPATCH_TABLE:
        if handler.software is not None:
            if not any(handler.software in str(s) for s in software):
                continue
        if handler.method is not None:
            if method != handler.method:
                continue
        if handler.algorithm is not None:
            if algorithm != handler.algorithm:
                continue
        if handler.dof is not None:
            if not handler.dof.issubset(dof):
                continue
        return handler
    return None


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


def _emit_input_params(ctx, step):
    """Emit stored input-parameter values as named variables and return
    a list of ``"name=var"`` strings ready to splice into a call.

    Each ASMO InputParameter stored in the graph becomes a Python variable
    initialised to its stored value so the user can easily override it::

        pressure = 0.0       # Pressure (PA)
        temperature = 300.0  # Temperature (K)
    """
    extra_kwargs = []
    for p in step.get("input_parameters") or []:
        raw_label = p.get("label") or "param"
        value = p.get("value")
        unit = p.get("unit") or ""
        var = _sanitize(raw_label)
        # Emit once (make_var deduplicates via _var_counts)
        var = ctx.make_var(var)
        comment = f"  # {raw_label} ({unit})" if unit else f"  # {raw_label}"
        ctx.code(f"{var} = {repr(value)}{comment}")
        extra_kwargs.append(f"{var}={var}")
    return extra_kwargs


def _handle_simulation(provenance, ctx, step):
    """Generate code for a Simulation step using _DISPATCH_TABLE."""
    in_id = step.get("input_sample_id")
    out_id = step["output_sample_id"]
    in_var = ctx.uri_to_var.get(in_id) if in_id else None

    if in_var is None:
        # Root step — no input, just load the output structure
        _ensure_structure(ctx, out_id, step.get("output_sample"))
        _register_calc_properties(provenance, ctx, step, var_map=None)
        return

    method = step.get("method") or "Simulation"
    out_label = _sample_label(ctx.kg, out_id)
    handler = _match_handler(step)

    if handler is None or handler.func is None:
        # No matching handler or handler has no implementation yet
        note = handler.note if handler else f"{method} (unrecognised software/method)"
        ctx.comment(note)
        ctx.code(f"# TODO: implement {note}")
        _ensure_structure(ctx, out_id, step.get("output_sample"))
        _register_calc_properties(provenance, ctx, step, var_map=None)
        return

    # Register user-fillable parameters (pair_style, pair_coeff, etc.)
    for param, desc in handler.user_inputs.items():
        ctx.add_user_input(param, desc)

    ctx.add_import(handler.import_line)

    # Emit stored input parameters as overridable variables
    extra_kwargs = _emit_input_params(ctx, step)

    # Build full kwargs string: handler defaults + graph-sourced params
    all_kwargs = handler.call_kwargs
    if extra_kwargs:
        all_kwargs = (
            all_kwargs + ", " + ", ".join(extra_kwargs)
            if all_kwargs
            else ", ".join(extra_kwargs)
        )

    ecoh_var = ctx.make_var("ecoh")
    vol_var = ctx.make_var("vol")

    if handler.returns_structure:
        out_var = ctx.make_var(f"atoms_{out_label}", uri=out_id)
        ctx.comment(f"{method} \u2192 {out_label} ({handler.note})")
        ctx.code(f"{out_var}, {ecoh_var}, {vol_var} = {handler.func}(")
        ctx.code(f"    {in_var}, {all_kwargs}")
        ctx.code(")")
    else:
        ctx.comment(f"{method} \u2192 {out_label} ({handler.note})")
        ctx.code(f"{ecoh_var}, {vol_var} = {handler.func}(")
        ctx.code(f"    {in_var}, {all_kwargs}")
        ctx.code(")")
        # Rigid: output sample is the same atoms object as input
        ctx.uri_to_var[out_id] = in_var

    _register_calc_properties(
        provenance,
        ctx,
        step,
        var_map={"energy": ecoh_var, "volume": vol_var},
    )


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
