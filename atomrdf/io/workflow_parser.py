import yaml
import json
import ast
from pathlib import Path
from typing import Dict, List, Union, Any, Optional
import re

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
from atomrdf import KnowledgeGraph
from rdflib import URIRef, Literal, Namespace, XSD

DCAT = Namespace("http://www.w3.org/ns/dcat#")

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

    Attributes
    ----------
    kg : KnowledgeGraph
        The knowledge graph to populate
    precision : int
        Decimal precision for hash computation
    sample_map : dict
        Maps original sample IDs to resolved URIs
    debug : bool
        If True, print debug messages during parsing
    hash_threshold : int
        Skip hashing for samples with more than this many atoms
    """

    def __init__(
        self,
        kg: Optional[KnowledgeGraph] = None,
        precision: int = 6,
        debug: bool = False,
        hash_threshold: int = 10000,
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
        hash_threshold : int, optional
            Skip hashing for samples with more than this many atoms.
            Default is 10000. Set to 0 to always hash.
        """
        self.kg = kg if kg is not None else KnowledgeGraph()
        self.precision = precision
        self.debug = debug
        self.hash_threshold = hash_threshold
        self.sample_map: Dict[str, str] = {}

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

    def parse_samples(self, sample_data_list: List[Dict[str, Any]]) -> Dict[str, str]:
        """
        Parse computational sample data and add to knowledge graph.

        Performs deduplication via hash-based lookup. If a sample with the same
        hash already exists, reuses the existing URI.

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
            sample_start = time.time()
            original_id = sample_data.get("id", "unknown")
            if self.debug:
                print(f"\n{'='*60}")
                print(f"Processing sample: {original_id}")

            # Normalise simulation_cell.vector if present
            prep_start = time.time()
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

            prep_time = time.time() - prep_start
            if self.debug:
                print(f"  Data preparation: {prep_time:.3f}s")

            # Create sample object
            model_start = time.time()
            sample = AtomicScaleSample(**sample_data)
            original_id = sample.id
            model_time = time.time() - model_start
            if self.debug:
                print(f"  Model creation: {model_time:.3f}s")

            # Check if we should skip hashing for large systems
            n_atoms = (
                simcell.get("number_of_atoms", 0) if isinstance(simcell, dict) else 0
            )
            skip_hash = self.hash_threshold > 0 and n_atoms > self.hash_threshold

            if skip_hash:
                # Skip hashing for large systems - treat as unique
                if self.debug:
                    print(
                        f"  Skipping hash for large system ({n_atoms} atoms > {self.hash_threshold} threshold)"
                    )
                graph_start = time.time()
                sample.id = None  # Let to_graph generate a new UUID
                sample.to_graph(self.kg)
                graph_time = time.time() - graph_start
                if self.debug:
                    print(f"  Graph addition: {graph_time:.3f}s")
                self.sample_map[original_id] = sample.id
            else:
                # Use hash-based deduplication for smaller systems
                hash_start = time.time()
                sample.id = None
                sample_hash = sample._compute_hash(precision=self.precision)
                hash_time = time.time() - hash_start
                if self.debug:
                    print(f"  Hash computation: {hash_time:.3f}s ({n_atoms} atoms)")
                    print(f"  Hash: {sample_hash}")

                # Check if this hash already exists in the KG
                lookup_start = time.time()
                existing_uri = self._find_sample_by_hash(sample_hash)
                lookup_time = time.time() - lookup_start
                if self.debug:
                    print(f"  Hash lookup: {lookup_time:.3f}s")
                    print(f"  Existing uri: {existing_uri}")

                if existing_uri:
                    self.sample_map[original_id] = existing_uri
                    if self.debug:
                        print(f"  Using existing sample (duplicate found)")
                else:
                    graph_start = time.time()
                    sample.to_graph(self.kg)
                    graph_time = time.time() - graph_start
                    if self.debug:
                        print(f"  Graph addition: {graph_time:.3f}s")
                    self.sample_map[original_id] = sample.id

            total_time = time.time() - sample_start
            if self.debug:
                print(f"  TOTAL sample time: {total_time:.3f}s")
                print(f"{'='*60}")

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

            # Create the Simulation object
            sim = Simulation(**workflow_data)

            # Add to knowledge graph
            sim_uri = sim.to_graph(self.kg)
            workflow_uris.append(sim_uri)

            if self.debug:
                print(
                    f"Added workflow {i+1}: connecting samples "
                    f"{workflow_data.get('input_sample', [])} to {workflow_data.get('output_sample', [])}"
                )

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
                print(
                    f"Added operation {i+1} ({method}): "
                    f"{operation_data.get('input_sample')} -> "
                    f"{operation_data.get('output_sample')}"
                )

        return operation_uris

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
            filepath = Path(data)

            if self.debug:
                print(f"\nReading file: {filepath}")

            read_start = time.time()
            with open(filepath, "r") as f:
                if filepath.suffix in [".yaml", ".yml"]:
                    data = yaml.safe_load(f)
                elif filepath.suffix == ".json":
                    data = json.load(f)
                else:
                    raise ValueError(f"Unsupported file format: {filepath.suffix}")
            read_time = time.time() - read_start

            if self.debug:
                print(f"File reading time: {read_time:.3f}s")
        elif not isinstance(data, dict):
            raise TypeError(
                f"Unsupported data type: {type(data)}. Expected str, Path, or dict."
            )

        result = {"sample_map": {}, "workflow_uris": [], "operation_uris": []}

        # Parse samples first (they may be referenced by workflows/operations)
        if "computational_sample" in data:
            result["sample_map"] = self.parse_samples(data["computational_sample"])

        # Parse workflows
        if "workflow" in data:
            result["workflow_uris"] = self.parse_workflows(data["workflow"])

        # Parse operations (formerly activities)
        if "operation" in data:
            result["operation_uris"] = self.parse_operations(data["operation"])
        # Backwards compatibility: also check for 'activity' key
        elif "activity" in data:
            result["operation_uris"] = self.parse_operations(data["activity"])

        return result

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
