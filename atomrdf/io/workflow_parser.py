import yaml
import json
import ast
from pathlib import Path
from typing import Dict, List, Union, Any, Optional
import re

from atomrdf.datamodels.dataset import Dataset
from atomrdf.datamodels.structure import AtomicScaleSample
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
    MATH_OPERATION_MAP,
)
from atomrdf import KnowledgeGraph
from rdflib import URIRef, Literal, Namespace, XSD

DCAT = Namespace("http://www.w3.org/ns/dcat#")

# Rebuild models to resolve forward references now that KnowledgeGraph is imported
DeleteAtom.model_rebuild()
SubstituteAtom.model_rebuild()
AddAtom.model_rebuild()
Rotate.model_rebuild()
Translate.model_rebuild()
Shear.model_rebuild()
Subtraction.model_rebuild()
Addition.model_rebuild()
Multiplication.model_rebuild()
Division.model_rebuild()
Exponentiation.model_rebuild()

# Mapping of operation method names to their classes
OPERATION_MAP = {
    "DeleteAtom": DeleteAtom,
    "SubstituteAtom": SubstituteAtom,
    "AddAtom": AddAtom,
    "Rotate": Rotate,
    "Rotation": Rotate,  # Alias
    "Translate": Translate,
    "Translation": Translate,  # Alias
    "Shear": Shear,
}


class WorkflowParser:
    """
    Parser for workflow YAML/JSON files into RDF knowledge graph.

    Handles parsing of:
    - Computational samples (with deduplication via hashing)
    - Workflows/Simulations
    - Operations (transformations between samples: DeleteAtom, SubstituteAtom,
      AddAtom, Rotate, Translate, Shear)
    - Math operations (ASMO arithmetic: Subtraction, Addition, Multiplication,
      Division, Exponentiation)

    Attributes
    ----------
    kg : KnowledgeGraph
        The knowledge graph to populate
    precision : int
        Decimal precision for hash computation
    sample_map : dict
        Maps original sample IDs to resolved URIs
    property_map : dict
        Maps user-defined property IDs (from YAML 'id' fields on
        calculated_property / input_parameter / output_parameter entries) to
        their generated KG URI strings.  Built incrementally as workflows and
        math operations are parsed, so later math_operation entries can
        reference earlier properties by their local ID.
    debug : bool
        If True, print debug messages during parsing
    hash_threshold : int or None
        Skip hashing for samples with more than this many atoms. Set to None to disable hashing completely.
    """

    def __init__(
        self,
        kg: Optional[KnowledgeGraph] = None,
        precision: int = 6,
        debug: bool = False,
        hash_threshold: Optional[int] = 10000,
    ):
        """
        Initialize the workflow parser.

        Parameters
        ----------
        kg : KnowledgeGraph, optional
            Knowledge graph instance. If None, creates a new one.
        precision : int, optional
            Decimal precision for sample hash computation. Default is 6.
        debug : bool, optional
            If True, print debug messages during parsing. Default is False.
        hash_threshold : int or None, optional
            Skip hashing for samples with more than this many atoms.
            Default is 10000. Set to 0 to always hash. Set to None to disable hashing completely.
        """
        self.kg = kg if kg is not None else KnowledgeGraph()
        self.precision = precision
        self.debug = debug
        self.hash_threshold = hash_threshold
        self.sample_map: Dict[str, str] = {}
        self.property_map: Dict[str, str] = {}
        self._base_dir: Optional[Path] = None
        # Raw YAML dicts keyed by sample id, used for dotpath resolution
        self._raw_sample_data: Dict[str, dict] = {}
        self._dotpath_cache: Dict[str, str] = {}

    def _miller_bravais_to_cartesian(self, uvtw: List[float]) -> List[float]:
        """
        Convert 4D Miller-Bravais indices [u, v, t, w] to 3D Cartesian [U, V, W].

        For hexagonal close-packed (hcp) systems, the Miller-Bravais notation
        uses four indices [u, v, t, w] where u + v + t = 0. This converts to
        3D Cartesian coordinates.

        Parameters
        ----------
        uvtw : list of float
            4-component Miller-Bravais indices [u, v, t, w]

        Returns
        -------
        list of float
            3-component Cartesian indices [U, V, W]

        Notes
        -----
        Conversion formula:
        U = 2u + v
        V = 2v + u
        W = w

        The third index t is redundant since t = -(u + v)
        """
        if len(uvtw) != 4:
            return uvtw

        u, v, t, w = uvtw
        # Convert to 3D using standard transformation
        U = 2 * u + v
        V = 2 * v + u
        W = w

        return [U, V, W]

    def _normalise_vector(self, vec: Any) -> Any:
        """
        Normalise a vector value into a 3x3 list-of-lists.

        Handles:
        - Already normalised 3x3 nested lists
        - 4D Miller-Bravais indices → 3D Cartesian conversion
        - Simple [x, y, z] lists → diagonal matrix
        - String representations of vectors

        Parameters
        ----------
        vec : Any
            Vector in various formats (list, nested list, or string)

        Returns
        -------
        Any
            Normalised 3x3 matrix or original value if cannot normalise
        """
        # Already a 3x3 matrix
        if (
            isinstance(vec, list)
            and len(vec) == 3
            and all(isinstance(v, list) for v in vec)
        ):
            # Check if any of the inner vectors are 4D (Miller-Bravais)
            converted_vec = []
            for v in vec:
                if len(v) == 4:
                    # Convert 4D Miller-Bravais to 3D
                    converted_vec.append(self._miller_bravais_to_cartesian(v))
                else:
                    converted_vec.append(v)
            return converted_vec

        # Simple [x, y, z] → diagonal matrix
        if (
            isinstance(vec, list)
            and len(vec) == 3
            and all(not isinstance(v, list) for v in vec)
        ):
            try:
                nums = [float(x) for x in vec]
                return [[nums[0], 0.0, 0.0], [0.0, nums[1], 0.0], [0.0, 0.0, nums[2]]]
            except Exception:
                return vec

        # String representation
        if isinstance(vec, str):
            # Handle ASE Cell object string representations
            # e.g., "Cell([[0.0, 1.785, 1.785], [1.785, 0.0, 1.785], [1.785, 1.785, 0.0]])"
            # or "Cell([2.862593415622699, 2.862593415622699, 2.862593415622699])"
            if vec.startswith("Cell("):
                # Extract content between Cell( and )
                inner = vec[5:-1]  # Remove "Cell(" and ")"
                try:
                    # Use ast.literal_eval to parse the list structure safely
                    parsed = ast.literal_eval(inner)

                    # Check if it's already a 3x3 matrix
                    if (
                        isinstance(parsed, list)
                        and len(parsed) == 3
                        and all(
                            isinstance(row, list) and len(row) == 3 for row in parsed
                        )
                    ):
                        return parsed

                    # Check if it's a simple [x, y, z] list
                    if (
                        isinstance(parsed, list)
                        and len(parsed) == 3
                        and all(not isinstance(v, list) for v in parsed)
                    ):
                        nums = [float(x) for x in parsed]
                        return [
                            [nums[0], 0.0, 0.0],
                            [0.0, nums[1], 0.0],
                            [0.0, 0.0, nums[2]],
                        ]
                except (ValueError, SyntaxError):
                    pass

            # Fallback: extract numbers from string
            nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", vec)
            if len(nums) == 3:
                numsf = [float(n) for n in nums]
                return [
                    [numsf[0], 0.0, 0.0],
                    [0.0, numsf[1], 0.0],
                    [0.0, 0.0, numsf[2]],
                ]
            elif len(nums) == 9:
                # Could be a flattened 3x3 matrix
                numsf = [float(n) for n in nums]
                return [
                    [numsf[0], numsf[1], numsf[2]],
                    [numsf[3], numsf[4], numsf[5]],
                    [numsf[6], numsf[7], numsf[8]],
                ]

        return vec

    def _convert_miller_bravais_in_dict(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Recursively convert 4D Miller-Bravais indices to 3D in dictionary fields.

        Specifically handles defect-related fields like:
        - stacking_fault.plane
        - grain_boundary fields
        - Any other vector-like fields

        Parameters
        ----------
        data : dict
            Dictionary potentially containing 4D vectors

        Returns
        -------
        dict
            Dictionary with 4D vectors converted to 3D
        """
        if not isinstance(data, dict):
            return data

        converted = {}
        for key, value in data.items():
            if isinstance(value, dict):
                # Recursively process nested dictionaries
                converted[key] = self._convert_miller_bravais_in_dict(value)
            elif isinstance(value, list) and len(value) == 4:
                # Check if this is a 4D vector (all numeric)
                try:
                    nums = [float(x) for x in value]
                    # Convert 4D Miller-Bravais to 3D
                    converted[key] = self._miller_bravais_to_cartesian(nums)
                except (ValueError, TypeError):
                    # Not a numeric vector, keep as is
                    converted[key] = value
            else:
                converted[key] = value

        return converted

    def _find_sample_by_hash(self, sample_hash: str) -> Optional[str]:
        """
        Find an existing sample in the knowledge graph by its hash.

        Parameters
        ----------
        sample_hash : str
            The computed hash of the sample

        Returns
        -------
        str or None
            URI of the existing sample, or None if not found
        """
        for g in self.kg.graph.triples(
            (None, DCAT.checksum, Literal(sample_hash, datatype=XSD.string))
        ):
            return g[0].toPython()
        return None

    def _resolve_atom_attribute_from_file(
        self, sample_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        If atom_attribute contains a file_path key, read the structure file via
        ASE and replace atom_attribute with inline position/species data.
        Also fills simulation_cell from the file if it is absent or empty in the
        YAML.

        Parameters
        ----------
        sample_data : dict
            A single sample dictionary (may be mutated by reference via copy).

        Returns
        -------
        dict
            Shallow-copied sample_data with atom_attribute resolved to inline
            arrays (and simulation_cell populated if needed).
        """
        from ase.io import read as ase_read

        aa = sample_data.get("atom_attribute")
        if not isinstance(aa, dict):
            return sample_data

        file_path = aa.get("file_path")
        if file_path is None:
            return sample_data

        # Resolve path relative to the YAML's parent directory
        base = self._base_dir if self._base_dir is not None else Path.cwd()
        resolved = (base / file_path).resolve()
        if not resolved.exists():
            raise FileNotFoundError(
                f"Structure file not found: {resolved} "
                f"(file_path='{file_path}', base='{base}')"
            )

        # Build ASE read kwargs
        read_kwargs: Dict[str, Any] = {}
        fmt = aa.get("file_format", "")
        if fmt:
            read_kwargs["format"] = fmt
        if aa.get("file_species"):
            species_list = list(aa["file_species"])
            if fmt == "lammps-data":
                # lammps-data doesn't accept specorder; use Z_of_type instead.
                # Read the number of atom types from the file header.
                import re
                from ase.data import atomic_numbers as _anz

                n_types = 1
                with open(str(resolved)) as _fh:
                    for _line in _fh:
                        _m = re.match(r"\s*(\d+)\s+atom\s+types", _line)
                        if _m:
                            n_types = int(_m.group(1))
                            break
                read_kwargs["Z_of_type"] = {
                    i: _anz[species_list[min(i - 1, len(species_list) - 1)]]
                    for i in range(1, n_types + 1)
                }
            else:
                read_kwargs["specorder"] = species_list

        atoms = ase_read(str(resolved), **read_kwargs)

        # Shallow-copy so we don't mutate the caller's dict
        sample_data = dict(sample_data)
        sample_data["atom_attribute"] = {
            "position": atoms.get_positions().tolist(),
            "species": atoms.get_chemical_symbols(),
        }

        # Fill simulation_cell from file only when absent / empty in the YAML
        simcell = sample_data.get("simulation_cell")
        needs_fill = not isinstance(simcell, dict) or not any(
            simcell.get(k) for k in ("length", "vector", "number_of_atoms")
        )
        if needs_fill:
            sample_data["simulation_cell"] = {
                "volume": {"value": float(atoms.get_volume())},
                "number_of_atoms": len(atoms),
                "length": atoms.cell.lengths().tolist(),
                "vector": atoms.get_cell().tolist(),
                "angle": atoms.cell.angles().tolist(),
            }

        return sample_data

    def parse_samples(self, sample_data_list: List[Dict[str, Any]]) -> Dict[str, str]:
        """
        Parse computational sample data and add to knowledge graph.

        Performs deduplication via hash-based lookup. If a sample with the same
        hash already exists, reuses the existing URI. If atom_attribute contains
        a file_path key, the structure file is read via ASE and atoms are
        resolved before building the sample object.

        Parameters
        ----------
        sample_data_list : list of dict
            List of sample dictionaries

        Returns
        -------
        dict
            Dictionary mapping original sample IDs to resolved URIs
        """
        import time

        for sample_data in sample_data_list:
            # Resolve file_path references in atom_attribute before any processing
            sample_data = self._resolve_atom_attribute_from_file(sample_data)
            original_id = sample_data.get("id", "unknown")
            # Store a copy of the raw YAML dict for later dotpath resolution
            self._raw_sample_data[original_id] = sample_data

            # Normalise simulation_cell.vector if present
            simcell = sample_data.get("simulation_cell")
            if isinstance(simcell, dict):
                if "vector" in simcell:
                    simcell["vector"] = self._normalise_vector(simcell["vector"])

                # Convert repetitions from tuple/string to list if needed
                if "repetitions" in simcell:
                    reps = simcell["repetitions"]
                    # If it's a string representation like "(15, 15, 15)", parse it
                    if isinstance(reps, str):
                        reps = ast.literal_eval(reps)
                    # Convert tuple to list
                    if isinstance(reps, tuple):
                        reps = list(reps)
                    simcell["repetitions"] = reps

                # Flatten grains: move grain_size and number_of_grains to simulation_cell level
                if "grains" in simcell:
                    grains = simcell["grains"]
                    if isinstance(grains, dict):
                        if "grain_size" in grains:
                            simcell["grain_size"] = grains["grain_size"]
                        if "number_of_grains" in grains:
                            simcell["number_of_grains"] = grains["number_of_grains"]
                        # Remove the nested grains field
                        del simcell["grains"]

                # Update the sample_data with modified simcell
                sample_data["simulation_cell"] = simcell

            # Convert 4D Miller-Bravais indices in defect fields
            defect_fields = [
                "stacking_fault",
                "grain_boundary",
                "tilt_grain_boundary",
                "twist_grain_boundary",
                "symmetric_tilt_grain_boundary",
                "mixed_grain_boundary",
                "vacancy",
                "interstitial",
                "substitutional",
                "dislocation",
            ]

            for field in defect_fields:
                if field in sample_data:
                    sample_data[field] = self._convert_miller_bravais_in_dict(
                        sample_data[field]
                    )

            # Handle defect_complex field if present
            if "defect_complex" in sample_data:
                defect_complex_data = sample_data["defect_complex"]
                # Convert Miller-Bravais in defect_complex if needed
                sample_data["defect_complex"] = self._convert_miller_bravais_in_dict(
                    defect_complex_data
                )

            # Create sample object
            sample = AtomicScaleSample(**sample_data)
            original_id = sample.id

            # Check if we should skip hashing for large systems or if hashing is disabled
            n_atoms = (
                simcell.get("number_of_atoms", 0) if isinstance(simcell, dict) else 0
            )
            skip_hash = (self.hash_threshold is None) or (
                self.hash_threshold > 0 and n_atoms > self.hash_threshold
            )

            if skip_hash:
                # Skip hashing for large systems - treat as unique
                sample.id = None  # Let to_graph generate a new UUID
                sample.to_graph(self.kg)
                self.sample_map[original_id] = sample.id
                if self.debug:
                    print(f"Sample added (no hash check): {sample.id}")
            else:
                # Use hash-based deduplication for smaller systems
                sample.id = None
                sample_hash = sample._compute_hash(precision=self.precision)

                # Check if this hash already exists in the KG
                existing_uri = self._find_sample_by_hash(sample_hash)

                if existing_uri:
                    self.sample_map[original_id] = existing_uri
                    if self.debug:
                        print(f"Sample exists: {existing_uri}")
                else:
                    sample.to_graph(self.kg)
                    self.sample_map[original_id] = sample.id
                    if self.debug:
                        print(f"Sample added: {sample.id}")

        return self.sample_map

    def parse_workflows(self, workflow_data_list: List[Dict[str, Any]]) -> List[str]:
        """
        Parse workflow/simulation data and add to knowledge graph.

        Resolves sample references using the sample_map.

        Parameters
        ----------
        workflow_data_list : list of dict
            List of workflow dictionaries

        Returns
        -------
        list of str
            List of workflow URIs created
        """
        workflow_uris = []

        for i, workflow_data in enumerate(workflow_data_list):
            # Fix algorithm name: UniaxialTension -> TensileTest
            if (
                "algorithm" in workflow_data
                and workflow_data["algorithm"] == "UniaxialTension"
            ):
                workflow_data["algorithm"] = "TensileTest"

            # Resolve input sample references
            if "input_sample" in workflow_data:
                for j, inp in enumerate(workflow_data["input_sample"]):
                    if inp in self.sample_map:
                        workflow_data["input_sample"][j] = self.sample_map[inp]

            # Resolve output sample references
            if "output_sample" in workflow_data:
                for j, outp in enumerate(workflow_data["output_sample"]):
                    if outp in self.sample_map:
                        workflow_data["output_sample"][j] = self.sample_map[outp]

            if "calculated_property" in workflow_data:
                for count, prop in enumerate(workflow_data["calculated_property"]):
                    if "associate_to_sample" in prop:
                        sample_list = []
                        for a_sample in prop["associate_to_sample"]:
                            if a_sample in self.sample_map:
                                sample_list.append(self.sample_map[a_sample])
                            else:
                                sample_list.append(a_sample)

                        workflow_data["calculated_property"][count][
                            "associate_to_sample"
                        ] = sample_list

            if "output_parameter" in workflow_data:
                for count, prop in enumerate(workflow_data["output_parameter"]):
                    if "associate_to_sample" in prop:
                        sample_list = [
                            self.sample_map.get(s, s)
                            for s in prop["associate_to_sample"]
                        ]
                        workflow_data["output_parameter"][count][
                            "associate_to_sample"
                        ] = sample_list

            # Capture user-provided 'id' fields on properties so they can be
            # registered in property_map after to_graph assigns KG URIs.
            user_prop_ids: Dict[tuple, str] = {}
            for prop_key in (
                "calculated_property",
                "input_parameter",
                "output_parameter",
            ):
                for idx, prop_data in enumerate(workflow_data.get(prop_key, []) or []):
                    orig_id = (
                        prop_data.get("id") if isinstance(prop_data, dict) else None
                    )
                    if orig_id:
                        user_prop_ids[(prop_key, idx)] = orig_id

            # Normalise input_parameter entries from nested-key format
            # {ParamName: {value, unit, label, ...}} to flat basename format
            # {basename: ParamName, value: ..., unit: ..., label: ...}
            if "input_parameter" in workflow_data:
                normalised = []
                for entry in workflow_data["input_parameter"] or []:
                    if (
                        isinstance(entry, dict)
                        and "basename" not in entry
                        and len(entry) == 1
                    ):
                        key, val = next(iter(entry.items()))
                        flat = {"basename": key}
                        if isinstance(val, dict):
                            if key == "KpointMesh" and "type" in val:
                                flat["basename"] = val["type"]
                            for field in ("value", "unit", "label"):
                                if field in val:
                                    flat[field] = val[field]
                            if "width" in val and "value" not in val:
                                flat["value"] = val["width"]
                        elif val is not None:
                            flat["label"] = str(val)
                        normalised.append(flat)
                    else:
                        normalised.append(entry)
                workflow_data["input_parameter"] = normalised

            # Create the Simulation object
            sim = Simulation(**workflow_data)

            # Add to knowledge graph
            sim_uri = sim.to_graph(self.kg)
            workflow_uris.append(sim_uri)

            # Register property IDs in property_map now that to_graph has
            # assigned KG URIs (Property.to_graph sets self.id in-place).
            for prop_key, prop_attr in [
                ("calculated_property", "calculated_property"),
                ("input_parameter", "input_parameter"),
                ("output_parameter", "output_parameter"),
            ]:
                prop_list = getattr(sim, prop_attr, None) or []
                for idx, prop in enumerate(prop_list):
                    key = (prop_key, idx)
                    if key in user_prop_ids and prop.id:
                        self.property_map[user_prop_ids[key]] = prop.id

            if self.debug:
                print(f"Workflow added: {sim_uri}")

        return workflow_uris

    def parse_operations(self, operation_data_list: List[Dict[str, Any]]) -> List[str]:
        """
        Parse operation data (transformations between samples).

        Operations include: DeleteAtom, SubstituteAtom, AddAtom, Rotate,
        Translate, and Shear.

        Parameters
        ----------
        operation_data_list : list of dict
            List of operation dictionaries. Each must have:
            - 'method': The operation type (e.g., 'DeleteAtom', 'Rotate')
            - 'input_sample': Sample ID or list of sample IDs
            - 'output_sample': Sample ID or list of sample IDs
            - Additional method-specific parameters (e.g., rotation_matrix for Rotate)

        Returns
        -------
        list of str
            List of operation URIs created

        Raises
        ------
        ValueError
            If operation method is not recognized
        """
        operation_uris = []

        for i, operation_data in enumerate(operation_data_list):
            method = operation_data.get("method")
            if not method:
                if self.debug:
                    print(f"Skipping operation {i+1}: no method specified")
                continue

            # Get the operation class
            operation_class = OPERATION_MAP.get(method)
            if not operation_class:
                raise ValueError(
                    f"Unknown operation method: {method}. "
                    f"Available methods: {list(OPERATION_MAP.keys())}"
                )

            # Resolve input sample references
            if "input_sample" in operation_data:
                input_sample = operation_data["input_sample"]
                if isinstance(input_sample, list):
                    operation_data["input_sample"] = [
                        self.sample_map.get(s, s) for s in input_sample
                    ]
                else:
                    operation_data["input_sample"] = self.sample_map.get(
                        input_sample, input_sample
                    )

            # Resolve output sample references
            if "output_sample" in operation_data:
                output_sample = operation_data["output_sample"]
                if isinstance(output_sample, list):
                    operation_data["output_sample"] = [
                        self.sample_map.get(s, s) for s in output_sample
                    ]
                else:
                    operation_data["output_sample"] = self.sample_map.get(
                        output_sample, output_sample
                    )

            # Remove 'method' from data as it's not part of the model
            operation_data_copy = operation_data.copy()
            operation_data_copy.pop("method", None)

            # Create the operation object
            operation = operation_class(**operation_data_copy)

            # Add to knowledge graph
            operation.to_graph(self.kg)
            operation_uris.append(operation.id)

            if self.debug:
                print(f"Operation added ({method}): {operation.id}")

        return operation_uris

    def _resolve_dotpath(self, operand_str: str) -> Optional[str]:
        """Resolve a dotpath operand like ``Fe_ref_relaxed.simulation_cell.number_of_atoms``.

        If the first segment is a known sample id, walks the raw YAML data along
        the remaining segments (supports ``field[index]`` for list access) and
        creates a ``property:`` node in the KG linked to that sample via
        ``ASMO.hasCalculatedProperty``.  The URI is cached so repeated
        references return the same node.

        Parameters
        ----------
        operand_str:
            A dotted path string whose first segment is a sample id.

        Returns
        -------
        str or None
            The ``property:`` URI string, or ``None`` if the path cannot be
            resolved.
        """
        import re as _re
        from rdflib import (
            URIRef as _URIRef,
            Literal as _Literal,
            XSD as _XSD,
            RDFS as _RDFS,
        )

        # Return cached result if already resolved
        if operand_str in self._dotpath_cache:
            return self._dotpath_cache[operand_str]

        parts = operand_str.split(".")
        if len(parts) < 2:
            return None

        sample_id = parts[0]
        raw = self._raw_sample_data.get(sample_id)
        if raw is None:
            return None

        # Walk the remaining path segments
        current = raw
        for segment in parts[1:]:
            # Support field[index] syntax
            m = _re.match(r"^(\w+)\[(\d+)\]$", segment)
            if m:
                key, idx = m.group(1), int(m.group(2))
                if not isinstance(current, dict) or key not in current:
                    return None
                seq = current[key]
                if not isinstance(seq, (list, tuple)) or idx >= len(seq):
                    return None
                current = seq[idx]
            else:
                if not isinstance(current, dict) or segment not in current:
                    return None
                current = current[segment]

        # current should now be a scalar
        if not isinstance(current, (int, float)):
            return None

        value = float(current)

        # Build a human-readable label from the last path segment
        last_seg = _re.match(r"^(\w+)(\[\d+\])?$", parts[-1])
        label = (
            last_seg.group(1).replace("_", " ").title().replace(" ", "")
            if last_seg
            else parts[-1]
        )
        # e.g. "number_of_atoms" -> "NumberOfAtoms"

        # Create a property node in the KG
        prop_id = f"property:{label.lower()}_{__import__('uuid').uuid4()}"
        from atomrdf.namespace import ASMO as _ASMO, PROV as _PROV

        prop_uri = self.kg.create_node(prop_id, _ASMO.CalculatedProperty, label=label)
        self.kg.add(
            (_URIRef(prop_id), _ASMO.hasValue, _Literal(value, datatype=_XSD.float))
        )
        # Store the dotpath as a comment so code-generation can recover it
        self.kg.add(
            (
                _URIRef(prop_id),
                _RDFS.comment,
                _Literal(operand_str, datatype=_XSD.string),
            )
        )

        # Link to its owning sample
        sample_uri = self.sample_map.get(sample_id)
        if sample_uri:
            self.kg.add(
                (_URIRef(sample_uri), _ASMO.hasCalculatedProperty, _URIRef(prop_id))
            )

        # Register in property_map and cache
        self.property_map[operand_str] = prop_id
        self._dotpath_cache[operand_str] = prop_id

        if self.debug:
            print(f"Dotpath resolved: {operand_str} -> {prop_id} (value={value})")

        return prop_id

    def parse_math_operations(
        self, math_op_data_list: List[Dict[str, Any]]
    ) -> List[str]:
        """
        Parse math-operation entries (ASMO arithmetic activities).

        Each entry must have a ``type`` key (one of ``Subtraction``,
        ``Addition``, ``Multiplication``, ``Division``, ``Exponentiation``).
        Operands may be local property-ID strings (resolved via
        ``self.property_map``) or numeric scalars.  If the result carries an
        ``id`` field it is registered in ``property_map`` so subsequent
        math_operation entries can use it as an operand.

        Parameters
        ----------
        math_op_data_list : list of dict

        Returns
        -------
        list of str
            List of math-operation activity-ID strings created.

        Raises
        ------
        ValueError
            If the ``type`` field is missing or unrecognised.
        """
        math_op_uris = []

        for i, mo_data in enumerate(math_op_data_list):
            op_type = mo_data.get("type")
            if not op_type:
                if self.debug:
                    print(f"Skipping math_operation {i + 1}: no type specified")
                continue

            op_class = MATH_OPERATION_MAP.get(op_type)
            if not op_class:
                raise ValueError(
                    f"Unknown math operation type: {op_type}. "
                    f"Available types: {list(MATH_OPERATION_MAP.keys())}"
                )

            mo_data = dict(mo_data)  # don't mutate the original
            mo_data.pop("type", None)

            # Capture the result's user-provided id before building the object
            result_data = mo_data.get("result")
            result_user_id = (
                result_data.get("id") if isinstance(result_data, dict) else None
            )

            # Resolve single-value operand string references
            for operand_key in (
                "minuend",
                "subtrahend",
                "dividend",
                "divisor",
                "base",
                "exponent",
            ):
                if operand_key in mo_data and isinstance(mo_data[operand_key], str):
                    val = mo_data[operand_key]
                    # Try property_map first, then dotpath resolution
                    if val in self.property_map:
                        mo_data[operand_key] = self.property_map[val]
                    elif "." in val:
                        resolved = self._resolve_dotpath(val)
                        if resolved is not None:
                            mo_data[operand_key] = resolved

            # Resolve list operands
            for operand_key in ("addend", "factor"):
                if operand_key in mo_data:
                    resolved = []
                    for v in mo_data[operand_key]:
                        if isinstance(v, str):
                            if v in self.property_map:
                                resolved.append(self.property_map[v])
                            elif "." in v:
                                dp = self._resolve_dotpath(v)
                                resolved.append(dp if dp is not None else v)
                            else:
                                resolved.append(v)
                        else:
                            resolved.append(v)
                    mo_data[operand_key] = resolved

            # Resolve associate_to_sample on the result
            if isinstance(result_data, dict) and "associate_to_sample" in result_data:
                mo_data["result"] = dict(result_data)
                mo_data["result"]["associate_to_sample"] = [
                    self.sample_map.get(s, s)
                    for s in result_data["associate_to_sample"]
                ]

            math_op = op_class(**mo_data)
            op_uri = math_op.to_graph(self.kg)
            math_op_uris.append(op_uri)

            # Register result in property_map using the user's id
            if result_user_id and math_op.result and math_op.result.id:
                self.property_map[result_user_id] = math_op.result.id

            if self.debug:
                print(f"Math operation added ({op_type}): {op_uri}")

        return math_op_uris

    def parse(self, data: Union[str, Path, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Parse complete workflow data structure.

        Parameters
        ----------
        data : str, Path, or dict
            Either a file path (str/Path) to a YAML/JSON file, or a dictionary
            containing computational_sample, workflow, and/or activity keys

        Returns
        -------
        dict
            Dictionary with the following keys:

            - 'sample_map' : dict mapping original IDs to URIs
            - 'workflow_uris' : list of created workflow URIs
            - 'operation_uris' : list of created operation URIs

        Raises
        ------
        ValueError
            If file format is not supported (must be .yaml, .yml, or .json)
        TypeError
            If data type is not supported
        """
        import time

        # If data is a file path, read it first
        if isinstance(data, (str, Path)):
            filepath = Path(data).resolve()
            self._base_dir = filepath.parent

            with open(filepath, "r") as f:
                if filepath.suffix in [".yaml", ".yml"]:
                    data = yaml.safe_load(f)
                elif filepath.suffix == ".json":
                    data = json.load(f)
                else:
                    raise ValueError(f"Unsupported file format: {filepath.suffix}")
        elif not isinstance(data, dict):
            raise TypeError(
                f"Unsupported data type: {type(data)}. Expected str, Path, or dict."
            )

        result = {
            "sample_map": {},
            "workflow_uris": [],
            "operation_uris": [],
            "math_operation_uris": [],
            "dataset_uri": None,
        }

        # Parse samples first (they may be referenced by workflows/operations)
        if "computational_sample" in data:
            result["sample_map"] = self.parse_samples(data["computational_sample"])

        # Parse dataset metadata (resolve sample references using sample_map)
        if "dataset" in data:
            result["dataset_uri"] = self._parse_dataset(data["dataset"])

        # Parse workflows
        if "workflow" in data:
            result["workflow_uris"] = self.parse_workflows(data["workflow"])

        # Parse operations (formerly activities)
        if "operation" in data:
            result["operation_uris"] = self.parse_operations(data["operation"])
        # Backwards compatibility: also check for 'activity' key
        elif "activity" in data:
            result["operation_uris"] = self.parse_operations(data["activity"])

        # Parse math operations (processed after workflows so property_map is
        # already populated with properties produced by simulations)
        if "math_operation" in data:
            result["math_operation_uris"] = self.parse_math_operations(
                data["math_operation"]
            )

        return result

    def _parse_dataset(self, dataset_data: Dict[str, Any]) -> str:
        """
        Parse dataset metadata and add to knowledge graph.

        Parameters
        ----------
        dataset_data : dict
            Dataset dictionary from the top-level YAML ``dataset`` key.
            Must follow the dataset_template schema: identifier, title,
            creators, publication, samples.

        Returns
        -------
        str
            URI string of the created dcat:Dataset node.
        """
        # Resolve sample ID references using the already-built sample_map
        dataset_data = dict(dataset_data)
        dataset_data["samples"] = [
            self.sample_map.get(s, s) for s in dataset_data.get("samples", [])
        ]
        dataset = Dataset.from_dict(dataset_data)
        dataset_uri = dataset.to_graph(self.kg)
        if self.debug:
            print(f"Dataset added: {dataset_uri}")
        return str(dataset_uri)

    def from_file(self, filepath: Union[str, Path]) -> Dict[str, Any]:
        """
        Parse workflow data from a YAML or JSON file.

        This is a convenience method that calls parse() with a file path.

        Parameters
        ----------
        filepath : str or Path
            Path to YAML or JSON file

        Returns
        -------
        dict
            Parse results dictionary

        Raises
        ------
        ValueError
            If file format is not supported (must be .yaml, .yml, or .json)
        """
        return self.parse(filepath)


# Backwards-compatible functional API
def parse_workflow_yaml(
    yaml_data: Dict[str, Any], kg: KnowledgeGraph, precision: int = 6
) -> Dict[str, str]:
    """
    Parse workflow YAML data (backwards-compatible function).

    Parameters
    ----------
    yaml_data : dict
        Dictionary containing workflow data
    kg : KnowledgeGraph
        Knowledge graph to populate
    precision : int, optional
        Decimal precision for hash computation. Default is 6.

    Returns
    -------
    dict
        Dictionary mapping original sample IDs to URIs
    """
    parser = WorkflowParser(kg=kg, precision=precision)
    result = parser.parse(yaml_data)
    return result["sample_map"]


def from_workflow_input(
    data: Union[str, Path, Dict[str, Any]],
    graph: Optional[KnowledgeGraph] = None,
    precision: int = 6,
) -> Dict[str, Any]:
    """
    Main entry point for parsing workflow data.

    Parameters
    ----------
    data : str, Path, or dict
        File path (str/Path) or dictionary containing workflow data
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, creates a new one.
    precision : int, optional
        Decimal precision for hash computation. Default is 6.

    Returns
    -------
    dict
        Dictionary containing parsed results

    Raises
    ------
    TypeError
        If data type is not supported (must be str, Path, or dict)
    """
    parser = WorkflowParser(kg=graph, precision=precision)

    if isinstance(data, (str, Path)):
        return parser.from_file(data)
    elif isinstance(data, dict):
        return parser.parse(data)
    else:
        raise TypeError(f"Unsupported data type: {type(data)}")


def parse_sample(
    sample_data: Dict[str, Any], graph: Optional[KnowledgeGraph] = None
) -> str:
    """
    Parse a single computational sample.

    Parameters
    ----------
    sample_data : dict
        Sample dictionary
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, creates a new one.

    Returns
    -------
    str or None
        URI of the created/found sample
    """
    parser = WorkflowParser(kg=graph)
    result = parser.parse_samples([sample_data])
    return list(result.values())[0] if result else None


def parse_workflow(
    workflow_data: Dict[str, Any], graph: Optional[KnowledgeGraph] = None
) -> str:
    """
    Parse a single workflow.

    Parameters
    ----------
    workflow_data : dict
        Workflow dictionary
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, creates a new one.

    Returns
    -------
    str or None
        URI of the created workflow
    """
    parser = WorkflowParser(kg=graph)
    uris = parser.parse_workflows([workflow_data])
    return uris[0] if uris else None


def parse_operation(
    operation_data: Dict[str, Any], graph: Optional[KnowledgeGraph] = None
) -> str:
    """
    Parse a single operation.

    Parameters
    ----------
    operation_data : dict
        Operation dictionary with 'method' field specifying the operation type
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, creates a new one.

    Returns
    -------
    str or None
        URI of the created operation
    """
    parser = WorkflowParser(kg=graph)
    uris = parser.parse_operations([operation_data])
    return uris[0] if uris else None


def parse_generic(
    model_class: type, data: Dict[str, Any], graph: Optional[KnowledgeGraph] = None
) -> Any:
    """
    Generic parser for any Pydantic model with .to_graph() method.

    Parameters
    ----------
    model_class : type
        Pydantic model class
    data : dict
        Dictionary to parse
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, model is created but not added to graph.

    Returns
    -------
    Any
        Instance of the model class
    """
    instance = model_class(**data)
    if graph is not None and hasattr(instance, "to_graph"):
        instance.to_graph(graph)
    return instance
