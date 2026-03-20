"""
Provenance tracing for atomRDF knowledge graphs.
"""

from rdflib import URIRef, RDF

from atomrdf.namespace import ASMO, PROV
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.datamodels.structure_io import sample_to_ase
from atomrdf.datamodels.workflow.workflow import Simulation
from atomrdf.datamodels.workflow.operations import (
    DeleteAtom,
    SubstituteAtom,
    AddAtom,
    Rotate,
    Translate,
    Shear,
)
from atomrdf.datamodels.workflow.math_operations import (
    Subtraction,
    Addition,
    Multiplication,
    Division,
    Exponentiation,
)

# ------------------------------------------------------------------ #
# Activity type → from_graph class mapping                             #
# ------------------------------------------------------------------ #

_OPERATION_MAP = {
    str(ASMO.DeleteAtom): DeleteAtom,
    str(ASMO.SubstituteAtom): SubstituteAtom,
    str(ASMO.AddAtom): AddAtom,
    str(ASMO.Rotation): Rotate,
    str(ASMO.Translation): Translate,
    str(ASMO.Shear): Shear,
    # Math operations
    str(ASMO.Subtraction): Subtraction,
    str(ASMO.Addition): Addition,
    str(ASMO.Multiplication): Multiplication,
    str(ASMO.Division): Division,
    str(ASMO.Exponentiation): Exponentiation,
}

_SIMULATION_TYPES = {
    str(ASMO.EnergyCalculation),
    str(ASMO.Simulation),
}

_MATH_OP_TYPES = {
    str(ASMO.Subtraction),
    str(ASMO.Addition),
    str(ASMO.Multiplication),
    str(ASMO.Division),
    str(ASMO.Exponentiation),
}

_MATH_OPERAND_PREDICATES = [
    ASMO.hasMinuend,
    ASMO.hasSubtrahend,
    ASMO.hasAddend,
    ASMO.hasFactor,
    ASMO.hasDividend,
    ASMO.hasDivisor,
    ASMO.hasBase,
    ASMO.hasExponent,
]


def _uri(s):
    """Ensure *s* is an ``rdflib.URIRef``."""
    return s if isinstance(s, URIRef) else URIRef(str(s))


def _short(uri):
    """Human-readable short name from a URI."""
    s = str(uri)
    if "#" in s:
        return s.rsplit("#", 1)[-1]
    if ":" in s and "/" not in s.split(":")[-1]:
        return s
    return s.rsplit("/", 1)[-1]


def _short_label(uri):
    """Display label: type prefix + first 4 chars of UUID part.

    E.g. ``sample:17e70b03-...`` → ``sample:17e7``
    ``property:totalenergy_48cb2441-...`` → ``totalenergy:48cb``
    """
    import re

    s = str(uri)
    # Strip URL-style namespace, keep short name
    if "/" in s:
        s = s.rsplit("/", 1)[-1]
    # For "prefix:rest" style, abbreviate UUID in rest
    if ":" in s:
        prefix, rest = s.split(":", 1)
        uuid_match = re.search(r"([0-9a-f]{4,})", rest, re.I)
        if uuid_match:
            return f"{prefix}:{rest[:rest.index(uuid_match.group(1)) + 4]}"
        return f"{prefix}:{rest[:8]}"
    # For underscore-separated names like "totalenergy_48cb2441_..."
    uuid_match = re.search(r"[_-]([0-9a-f]{4,})", s, re.I)
    if uuid_match:
        return s[: s.index(uuid_match.group(1)) + 4]
    return s[:20]


class Provenance:
    """Trace the provenance of a sample or calculated property.

    Walks backwards through ``PROV.wasDerivedFrom`` /
    ``PROV.wasGeneratedBy`` links and uses the ``from_graph``
    methods of each data-model class to reconstruct every object.

    Parameters
    ----------
    kg : KnowledgeGraph
        The knowledge graph to query.

    Attributes
    ----------
    pipeline : list of dict
        Ordered list of step dicts
    """

    def __init__(self, kg):
        self.kg = kg
        self._steps = []
        self._sample_cache = {}  # uri-str → (AtomicScaleSample, Atoms|None)
        self._property_uri = None

    @classmethod
    def from_sample(cls, kg, sample_uri):
        """Trace provenance backwards from *sample_uri*.

        Parameters
        ----------
        kg : KnowledgeGraph
        sample_uri : str or URIRef

        Returns
        -------
        Provenance
        """
        prov = cls(kg)
        prov._trace(_uri(sample_uri))
        return prov

    @classmethod
    def from_property(cls, kg, property_uri):
        """Trace provenance backwards from a calculated property.

        Finds the sample that owns the property, then traces from there.

        Parameters
        ----------
        kg : KnowledgeGraph
        property_uri : str or URIRef

        Returns
        -------
        Provenance

        Raises
        ------
        ValueError
            If no sample in the graph owns the property.
        """
        prop = _uri(property_uri)
        prov = cls(kg)
        prov._property_uri = prop

        sample = None
        for s, _, _ in kg.triples((None, ASMO.hasCalculatedProperty, prop)):
            sample = s
            break
        if sample is None:
            raise ValueError(f"No sample found with property {property_uri}")

        prov._trace(sample)

        # Append math-operation steps that produced this property
        for act_uri, result_uri in prov._collect_math_ops(kg, prop):
            prov._steps.append(prov._build_math_step(kg, act_uri, result_uri))

        return prov

    @property
    def pipeline(self):
        """Ordered list of step dicts"""
        return list(self._steps)

    def __len__(self):
        return len(self._steps)

    def __iter__(self):
        return iter(self._steps)

    def __getitem__(self, idx):
        return self._steps[idx]

    def __repr__(self):
        types = " → ".join(s["activity_type"] for s in self._steps)
        return f"Provenance({len(self._steps)} steps: {types})"

    def visualize(self, rankdir="LR", size=None, layout="dot", filename=None, dpi=150):
        """Return a ``graphviz.Digraph`` of the provenance chain.

        Samples are shown as ellipses, activities / operations as
        boxes, and calculated properties as diamonds.

        Parameters
        ----------
        rankdir : str
            Layout direction (``"LR"``, ``"TB"``, …).
        size : tuple of int, optional
            ``(width, height)`` in inches.
        layout : str
            Graphviz layout engine.
        filename : str, optional
            If given, render and save the diagram to this path.  The file
            format is inferred from the extension (e.g. ``"prov.pdf"``,
            ``"prov.png"``, ``"prov.svg"``).  Defaults to PNG when no
            extension is recognised.
        dpi : int
            Resolution in dots-per-inch used when saving raster formats
            (PNG, …).  Ignored for vector formats (PDF, SVG).  Default 150.

        Returns
        -------
        graphviz.Digraph
        """
        import graphviz
        import re as _re

        from atomrdf.namespace import CMSO as _CMSO

        # Graphviz interprets "name:port" on colons in node IDs, so sanitize
        # all IDs used as graph node names.
        def _gvid(uri):
            return _re.sub(r"[^A-Za-z0-9_]", "_", str(uri))

        def _composition_label(sample_uri):
            """Build a label like 'Fe\n(bcc)' or 'Fe Al₀.₀₁\n(Im-3m)' from the KG."""
            species_node = self.kg.graph.value(_uri(sample_uri), _CMSO.hasSpecies)
            composition = {}
            if species_node is not None:
                for el in self.kg.graph.objects(species_node, _CMSO.hasElement):
                    sym = self.kg.graph.value(el, _CMSO.hasChemicalSymbol)
                    ratio = self.kg.graph.value(el, _CMSO.hasElementRatio)
                    if sym is not None:
                        composition[str(sym)] = (
                            float(ratio) if ratio is not None else 1.0
                        )

            if not composition:
                return self.kg.get_label(_uri(sample_uri)) or _short_label(
                    str(sample_uri)
                )

            # Sort by descending ratio so the host element comes first
            sorted_els = sorted(composition.items(), key=lambda x: -x[1])
            if len(sorted_els) == 1:
                comp_str = sorted_els[0][0]
            else:
                parts = []
                for el, r in sorted_els:
                    if r > 0.999:
                        parts.append(el)
                    else:
                        parts.append(f"{el}{r:.3f}")
                comp_str = " ".join(parts)

            return comp_str

        dot = graphviz.Digraph()
        dot.attr(rankdir=rankdir, layout=layout, overlap="false")
        if size is not None:
            dot.attr(size=f"{size[0]},{size[1]}")

        seen_nodes = set()
        seen_edges = set()

        def _add_edge(src, dst, **attrs):
            key = (src, dst)
            if key not in seen_edges:
                seen_edges.add(key)
                dot.edge(src, dst, **attrs)

        def _ensure_sample_node(sid):
            if sid and sid not in seen_nodes:
                seen_nodes.add(sid)
                slabel = _composition_label(sid)
                dot.node(
                    _gvid(sid),
                    label=slabel,
                    shape="ellipse",
                    style="filled",
                    color="#C9DAF8",
                    fontname="Helvetica",
                    fontsize="9",
                )

        def _render_chain_steps(steps_list):
            """Render steps as sample nodes with activity labels on edges."""
            for s in steps_list:
                in_id = s.get("input_sample_id")
                out_id = s.get("output_sample_id")
                _ensure_sample_node(in_id)
                _ensure_sample_node(out_id)

                # Build edge label from method / activity type
                edge_label = s.get("method") or s["activity_type"]
                if s.get("algorithm"):
                    edge_label += f"\n({s['algorithm']})"

                if in_id and out_id:
                    _add_edge(
                        _gvid(in_id),
                        _gvid(out_id),
                        label=edge_label,
                        color="#263238",
                        fontname="Helvetica",
                        fontsize="8",
                        fontcolor="#B22222",
                    )
                elif out_id:
                    # Root step with no input — just ensure the node exists
                    pass

        _render_chain_steps([s for s in self._steps if "result_property" not in s])

        # Math operation steps — rendered as property → activity → property
        math_steps = [s for s in self._steps if "result_property" in s]

        for i, step in enumerate(math_steps):
            act_id = step["activity_id"]
            act_gv = _gvid(act_id)
            rp = step["result_property"]
            rp_id = rp["uri"]
            rp_gv = _gvid(rp_id)
            rp_label = rp["label"] or _short_label(rp_id)
            if rp["value"] is not None:
                rp_label += f"\n= {rp['value']:.4g}"
                if rp["unit"]:
                    rp_label += f" {rp['unit']}"

            # Activity node (yellow box) — kept for math ops only
            if act_id not in seen_nodes:
                seen_nodes.add(act_id)
                dot.node(
                    act_gv,
                    label=step["activity_type"],
                    shape="box",
                    style="filled",
                    color="#FFE599",
                    fontname="Helvetica",
                    fontsize="9",
                )
            # Result property diamond
            if rp_id not in seen_nodes:
                seen_nodes.add(rp_id)
                dot.node(
                    rp_gv,
                    label=rp_label,
                    shape="diamond",
                    style="filled",
                    color="#D5E8D4",
                    fontname="Helvetica",
                    fontsize="8",
                )
            _add_edge(
                act_gv, rp_gv, color="#263238", fontname="Helvetica", fontsize="7"
            )

            # Connect operand properties to this activity
            if step["activity"] is not None:
                op = step["activity"]
                operand_attrs = []
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
                    val = getattr(op, attr, None)
                    if val is None:
                        continue
                    items = val if isinstance(val, list) else [val]
                    operand_attrs.extend(items)
                for item in operand_attrs:
                    if isinstance(item, str) and item.startswith("property:"):
                        # Ensure the operand node exists in the graph
                        if item not in seen_nodes:
                            seen_nodes.add(item)
                            # Use rdfs:comment (dotpath) as label if present
                            from rdflib import RDFS as _RDFS

                            comment_n = self.kg.graph.value(
                                _uri(item),
                                _RDFS.comment,
                            )
                            if comment_n is not None:
                                # e.g. "Fe_ref_relaxed.simulation_cell.number_of_atoms"
                                # show only the last segment(s) for brevity
                                dotpath = str(comment_n)
                                parts = dotpath.split(".")
                                op_label = (
                                    ".".join(parts[1:]) if len(parts) > 1 else dotpath
                                )
                            else:
                                op_label = self.kg.get_label(
                                    _uri(item)
                                ) or _short_label(item)
                            # Try to get value/unit for display
                            val_n = self.kg.graph.value(_uri(item), ASMO.hasValue)
                            unit_n = self.kg.graph.value(_uri(item), ASMO.hasUnit)
                            if val_n is not None:
                                op_label += f"\n= {float(val_n):.4g}"
                                if unit_n:
                                    op_label += f" {str(unit_n).rsplit('/', 1)[-1]}"
                            dot.node(
                                _gvid(item),
                                label=op_label,
                                shape="diamond",
                                style="filled",
                                color="#D5E8D4",
                                fontname="Helvetica",
                                fontsize="8",
                            )
                            # If this property is owned by a sample, draw a
                            # dotted edge from that sample to the property
                            # diamond.  If the owner is not yet in the graph
                            # (out-of-chain reference sample), trace its
                            # provenance and render those steps so the
                            # simulation that produced it is also visible.
                            # If the sub-prov is empty (bare root sample with
                            # no simulation history), skip creating a new
                            # disconnected node — a bare orphan would cause a
                            # spurious second top-level node in the diagram.
                            owner = self.kg.graph.value(
                                None, ASMO.hasCalculatedProperty, _uri(item)
                            )
                            if owner is not None:
                                owner_str = str(owner)
                                if owner_str not in seen_nodes:
                                    # Render provenance of this reference sample
                                    _owner_prov = Provenance.from_sample(
                                        self.kg, owner
                                    )
                                    _render_chain_steps(
                                        [
                                            s
                                            for s in _owner_prov._steps
                                            if "result_property" not in s
                                        ]
                                    )
                                    # Only create a bare node as fallback when
                                    # the sub-prov produced at least one step
                                    # (i.e. the owner has some simulation
                                    # history).  If _steps is empty the owner
                                    # is a pristine imported structure with no
                                    # provenance — adding it as a bare orphan
                                    # would create a spurious top-level node.
                                    if _owner_prov._steps:
                                        _ensure_sample_node(owner_str)
                                if owner_str in seen_nodes:
                                    _add_edge(
                                        _gvid(owner_str),
                                        _gvid(item),
                                        color="#888888",
                                        fontname="Helvetica",
                                        fontsize="7",
                                        style="dotted",
                                    )
                        _add_edge(
                            _gvid(item),
                            act_gv,
                            color="#263238",
                            fontname="Helvetica",
                            fontsize="7",
                            style="dashed",
                        )
                    elif isinstance(item, (int, float)):
                        # Numeric literal — create a small constant node
                        const_id = f"const_{act_id}_{item}"
                        const_gv = _gvid(const_id)
                        if const_id not in seen_nodes:
                            seen_nodes.add(const_id)
                            dot.node(
                                const_gv,
                                label=f"{item:g}",
                                shape="plaintext",
                                fontname="Helvetica",
                                fontsize="8",
                            )
                        _add_edge(
                            const_gv,
                            act_gv,
                            color="#263238",
                            fontname="Helvetica",
                            fontsize="7",
                            style="dashed",
                        )

        # If tracing from a property with no math steps, draw old-style diamond
        if self._property_uri:
            pid = str(self._property_uri)
            if not math_steps:
                plabel = self.kg.get_label(_uri(pid)) or _short_label(pid)
                val_n = self.kg.graph.value(_uri(pid), ASMO.hasValue)
                unit_n = self.kg.graph.value(_uri(pid), ASMO.hasUnit)
                if val_n is not None:
                    try:
                        plabel += f"\n= {float(val_n):.4g}"
                    except (ValueError, TypeError):
                        plabel += f"\n= {val_n}"
                    if unit_n:
                        plabel += f" {str(unit_n).rsplit('/', 1)[-1]}"
                if pid not in seen_nodes:
                    dot.node(
                        _gvid(pid),
                        label=plabel,
                        shape="diamond",
                        style="filled",
                        color="#D5E8D4",
                        fontname="Helvetica",
                        fontsize="8",
                    )
                last_sample_id = next(
                    (
                        s["output_sample_id"]
                        for s in reversed(self._steps)
                        if s.get("output_sample_id")
                    ),
                    None,
                )
                if last_sample_id:
                    _add_edge(
                        _gvid(last_sample_id),
                        _gvid(pid),
                        label="hasCalculatedProperty",
                        color="#263238",
                        fontname="Helvetica",
                        fontsize="7",
                    )

        if filename is not None:
            import os

            root, ext = os.path.splitext(filename)
            fmt = ext.lstrip(".").lower() if ext else "png"
            if fmt not in ("pdf", "svg", "png", "jpg", "jpeg", "eps"):
                fmt = "png"
            # DPI is passed as a graph attribute for raster formats
            if fmt not in ("pdf", "svg", "eps"):
                dot.attr(dpi=str(dpi))
            dot.render(filename=root, format=fmt, cleanup=True)

        return dot

    def to_code(self, output_dir=None):
        """Generate a Python workflow script from this provenance chain.

        Walks the traced steps and emits executable Python that uses
        real workflow functions (e.g. ``atomrdf.workflow.lammps``) to
        reproduce the computation.  Atom operations are not re-executed;
        their output structures are written as files instead.

        Parameters
        ----------
        output_dir : str, optional
            If given, write ``workflow.py`` and ``structures/`` here.

        Returns
        -------
        str
            The generated Python script.
        """
        from atomrdf.io.codegen import generate_code

        return generate_code(self, output_dir)

    # -- internals ---------------------------------------------------- #

    def _trace(self, sample_uri):
        """Walk backwards from *sample_uri* and build steps root-first."""
        raw = []  # [(input_uri | None, activity_uri, output_uri)]
        current = sample_uri
        visited = set()

        while current is not None and str(current) not in visited:
            visited.add(str(current))
            activity = self.kg.value(current, PROV.wasGeneratedBy)
            if activity is None:
                break
            parent = self.kg.value(current, PROV.wasDerivedFrom)
            raw.append((parent, activity, current))
            current = parent

        # Reverse so first step is the root (oldest ancestor)
        raw.reverse()

        for in_uri, act_uri, out_uri in raw:
            self._steps.append(self._build_step(in_uri, act_uri, out_uri))

    def _load_sample(self, uri):
        """Load and cache ``AtomicScaleSample`` + ASE ``Atoms``."""
        key = str(uri)
        if key in self._sample_cache:
            return self._sample_cache[key]

        sample_obj, atoms = None, None
        try:
            sample_obj = AtomicScaleSample.from_graph(self.kg, key)
            atoms = sample_to_ase(sample_obj)
        except Exception:
            pass

        self._sample_cache[key] = (sample_obj, atoms)
        return sample_obj, atoms

    def _classify_activity(self, act_uri):
        """Return ``(type_name, is_simulation, operation_class_or_None)``."""
        rdf_type = self.kg.value(_uri(act_uri), RDF.type)
        if rdf_type is None:
            return "Unknown", False, None
        type_str = str(rdf_type)
        if type_str in _SIMULATION_TYPES:
            return "Simulation", True, None
        if type_str in _OPERATION_MAP:
            op_cls = _OPERATION_MAP[type_str]
            return op_cls.__name__, False, op_cls
        return _short(rdf_type), False, None

    def _build_step(self, in_uri, act_uri, out_uri):
        """Construct a single pipeline step dict using ``from_graph``."""
        act_type, is_sim, op_cls = self._classify_activity(act_uri)

        # Load structures via AtomicScaleSample.from_graph
        if in_uri is not None:
            _, in_atoms = self._load_sample(in_uri)
        else:
            in_atoms = None
        _, out_atoms = self._load_sample(out_uri)

        step = {
            "activity_type": act_type,
            "activity_id": str(act_uri),
            "activity": None,
            "input_sample": in_atoms,
            "output_sample": out_atoms,
            "input_sample_id": str(in_uri) if in_uri else None,
            "output_sample_id": str(out_uri),
            "method": None,
            "algorithm": None,
            "input_parameters": [],
            "output_parameters": [],
            "calculated_properties": [],
            "interatomic_potential": None,
            "degrees_of_freedom": [],
            "software": [],
        }

        if is_sim:
            self._fill_from_simulation(self.kg, act_uri, step)
        elif op_cls is not None:
            self._fill_from_operation(self.kg, act_uri, op_cls, step)

        return step

    @staticmethod
    def _fill_from_simulation(kg, sim_uri, step):
        """Use ``Simulation.from_graph`` to populate the step dict."""
        sim = Simulation.from_graph(kg, str(sim_uri))
        step["activity"] = sim

        # Method
        method = getattr(sim, "method", None)
        if method is not None:
            step["method"] = (
                method.basename if hasattr(method, "basename") else str(method)
            )

        # Algorithm — prefer explicit algorithm, fall back to thermodynamic ensemble
        # (e.g. MolecularDynamics stores IsothermalIsobaricEnsemble under
        #  hasStatisticalEnsemble, which the codegen dispatch table matches as algorithm)
        algo = getattr(sim, "algorithm", None)
        if algo is not None:
            step["algorithm"] = (
                algo.basename if hasattr(algo, "basename") else str(algo)
            )
        else:
            ensemble = getattr(sim, "thermodynamic_ensemble", None)
            if ensemble is not None:
                step["algorithm"] = (
                    ensemble.basename
                    if hasattr(ensemble, "basename")
                    else str(ensemble).split("/")[-1]
                )

        # Input parameters
        for p in getattr(sim, "input_parameter", None) or []:
            step["input_parameters"].append(
                {
                    "label": getattr(p, "label", None),
                    "value": getattr(p, "value", None),
                    "unit": getattr(p, "unit", None),
                }
            )

        # Output parameters
        for p in getattr(sim, "output_parameter", None) or []:
            step["output_parameters"].append(
                {
                    "label": getattr(p, "label", None),
                    "value": getattr(p, "value", None),
                    "unit": getattr(p, "unit", None),
                }
            )

        # Calculated properties
        for p in getattr(sim, "calculated_property", None) or []:
            step["calculated_properties"].append(
                {
                    "label": getattr(p, "label", None),
                    "value": getattr(p, "value", None),
                    "unit": getattr(p, "unit", None),
                }
            )

        # Interatomic potential
        pot = getattr(sim, "interatomic_potential", None)
        if pot is not None and hasattr(pot, "potential_type"):
            step["interatomic_potential"] = {
                "type": pot.potential_type,
                "uri": getattr(pot, "uri", None),
                "label": getattr(pot, "label", None),
            }

        # Degrees of freedom
        for dof in getattr(sim, "degrees_of_freedom", None) or []:
            name = dof.basename if hasattr(dof, "basename") else str(dof)
            step["degrees_of_freedom"].append(name)

        # Software
        for sw in getattr(sim, "software", None) or []:
            name = getattr(sw, "label", None) or getattr(sw, "uri", None) or str(sw)
            step["software"].append(name)

    @staticmethod
    def _collect_math_ops(kg, prop_uri):
        """Walk backwards from *prop_uri* through math-op chains.

        Returns a list of ``(activity_uri, result_property_uri)`` pairs in
        **execution order** (leaves first, final op last) via DFS post-order.
        """
        from rdflib import Literal as RDFLiteral

        result = []
        visited = set()

        def dfs(prop):
            key = str(prop)
            if key in visited:
                return
            visited.add(key)
            activity = kg.graph.value(_uri(prop), PROV.wasGeneratedBy)
            if activity is None:
                return
            rdf_type = kg.graph.value(activity, RDF.type)
            if rdf_type is None or str(rdf_type) not in _MATH_OP_TYPES:
                return
            # Recurse into operand properties
            for pred in _MATH_OPERAND_PREDICATES:
                for operand in kg.graph.objects(activity, pred):
                    if isinstance(operand, URIRef) and str(operand).startswith(
                        "property:"
                    ):
                        dfs(operand)
            result.append((activity, prop))

        dfs(_uri(prop_uri))
        return result  # DFS post-order = execution order

    def _build_math_step(self, kg, act_uri, result_uri):
        """Build a pipeline step dict for a single math operation."""
        rdf_type = kg.graph.value(act_uri, RDF.type)
        act_type = _short(rdf_type) if rdf_type else "MathOperation"
        op_cls = _OPERATION_MAP.get(str(rdf_type))

        label_node = kg.graph.value(
            result_uri, URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        )
        val_node = kg.graph.value(result_uri, ASMO.hasValue)
        unit_node = kg.graph.value(result_uri, ASMO.hasUnit)

        result_prop = {
            "uri": str(result_uri),
            "label": str(label_node) if label_node else act_type,
            "value": float(val_node) if val_node is not None else None,
            "unit": str(unit_node).split("/")[-1] if unit_node else None,
        }

        step = {
            "activity_type": act_type,
            "activity_id": str(act_uri),
            "activity": op_cls.from_graph(kg, str(act_uri)) if op_cls else None,
            "input_sample": None,
            "output_sample": None,
            "input_sample_id": None,
            "output_sample_id": None,
            "result_property": result_prop,
            "method": None,
            "algorithm": None,
            "input_parameters": [],
            "output_parameters": [],
            "calculated_properties": [],
            "interatomic_potential": None,
            "degrees_of_freedom": [],
            "software": [],
        }
        return step

    @staticmethod
    def _fill_from_operation(kg, act_uri, op_cls, step):
        """Use an operation's ``from_graph`` to populate the step dict."""
        op = op_cls.from_graph(kg, str(act_uri))
        step["activity"] = op
