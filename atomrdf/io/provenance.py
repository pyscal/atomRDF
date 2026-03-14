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

    def visualize(self, rankdir="LR", size=None, layout="dot"):
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

        Returns
        -------
        graphviz.Digraph
        """
        import graphviz

        dot = graphviz.Digraph()
        dot.attr(rankdir=rankdir, layout=layout, overlap="false")
        if size is not None:
            dot.attr(size=f"{size[0]},{size[1]}")

        seen_nodes = set()

        for step in self._steps:
            # Activity node
            act_id = step["activity_id"]
            label = step["activity_type"]
            if step.get("method"):
                label = step["method"]
                if step.get("algorithm"):
                    label += f"\n({step['algorithm']})"

            dot.node(
                act_id,
                label=label,
                shape="box",
                style="filled",
                color="#E6B8AF",
                fontname="Helvetica",
                fontsize="9",
            )

            # Sample nodes
            for sid in (step.get("input_sample_id"), step["output_sample_id"]):
                if sid and sid not in seen_nodes:
                    seen_nodes.add(sid)
                    slabel = self.kg.get_label(_uri(sid)) or _short(sid)
                    dot.node(
                        sid,
                        label=slabel,
                        shape="ellipse",
                        style="filled",
                        color="#C9DAF8",
                        fontname="Helvetica",
                        fontsize="8",
                    )

            # Edges
            if step.get("input_sample_id"):
                dot.edge(
                    step["input_sample_id"],
                    act_id,
                    color="#263238",
                    fontname="Helvetica",
                    fontsize="7",
                )
            dot.edge(
                act_id,
                step["output_sample_id"],
                color="#263238",
                fontname="Helvetica",
                fontsize="7",
            )

        # Property node (if tracing from a calculated property)
        if self._property_uri and self._steps:
            pid = str(self._property_uri)
            plabel = self.kg.get_label(_uri(pid)) or _short(pid)
            dot.node(
                pid,
                label=plabel,
                shape="diamond",
                style="filled",
                color="#D5E8D4",
                fontname="Helvetica",
                fontsize="8",
            )
            dot.edge(
                self._steps[-1]["output_sample_id"],
                pid,
                label="hasCalculatedProperty",
                color="#263238",
                fontname="Helvetica",
                fontsize="7",
            )

        return dot

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

        # Algorithm
        algo = getattr(sim, "algorithm", None)
        if algo is not None:
            step["algorithm"] = (
                algo.basename if hasattr(algo, "basename") else str(algo)
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
            name = getattr(sw, "name", None) or str(sw)
            step["software"].append(name)

    @staticmethod
    def _fill_from_operation(kg, act_uri, op_cls, step):
        """Use an operation's ``from_graph`` to populate the step dict."""
        op = op_cls.from_graph(kg, str(act_uri))
        step["activity"] = op
