import yaml
import json
from pathlib import Path
from typing import Dict, List, Union, Any, Optional
import re

from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.datamodels.workflow.workflow import Simulation
from atomrdf import KnowledgeGraph
from rdflib import URIRef, Literal, Namespace, XSD

DCAT = Namespace("http://www.w3.org/ns/dcat#")


class WorkflowParser:
    """
    Parser for workflow YAML/JSON files into RDF knowledge graph.

    Handles parsing of:
    - Computational samples (with deduplication via hashing)
    - Workflows/Simulations
    - Activities (transformations between samples)

    Attributes
    ----------
    kg : KnowledgeGraph
        The knowledge graph to populate
    precision : int
        Decimal precision for hash computation
    sample_map : dict
        Maps original sample IDs to resolved URIs
    """

    def __init__(self, kg: Optional[KnowledgeGraph] = None, precision: int = 6):
        """
        Initialize the workflow parser.

        Parameters
        ----------
        kg : KnowledgeGraph, optional
            Knowledge graph instance. If None, creates a new one.
        precision : int, optional
            Decimal precision for sample hash computation. Default is 6.
        """
        self.kg = kg if kg is not None else KnowledgeGraph()
        self.precision = precision
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
            nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", vec)
            if len(nums) == 3:
                numsf = [float(n) for n in nums]
                return [
                    [numsf[0], 0.0, 0.0],
                    [0.0, numsf[1], 0.0],
                    [0.0, 0.0, numsf[2]],
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
        for sample_data in sample_data_list:
            original_id = sample_data.get("id", "unknown")
            print(f"Processing sample: {original_id}")

            # Normalise simulation_cell.vector if present
            simcell = sample_data.get("simulation_cell")
            if isinstance(simcell, dict) and "vector" in simcell:
                simcell["vector"] = self._normalise_vector(simcell["vector"])
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

            # Create sample object to compute hash
            sample = AtomicScaleSample(**sample_data)
            original_id = sample.id
            sample.id = None
            sample_hash = sample._compute_hash(precision=self.precision)
            print(f"Computed hash: {sample_hash}")

            # Check if this hash already exists in the KG
            existing_uri = self._find_sample_by_hash(sample_hash)
            print(f"Existing uri is {existing_uri}")

            if existing_uri:
                self.sample_map[original_id] = existing_uri
            else:
                sample.to_graph(self.kg)
                self.sample_map[original_id] = sample.id

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
            # Resolve input sample references
            if "inputs" in workflow_data:
                for j, inp in enumerate(workflow_data["inputs"]):
                    if inp in self.sample_map:
                        workflow_data["inputs"][j] = self.sample_map[inp]

            # Resolve output sample references
            if "outputs" in workflow_data:
                for j, outp in enumerate(workflow_data["outputs"]):
                    if outp in self.sample_map:
                        workflow_data["outputs"][j] = self.sample_map[outp]

            # Create the Simulation object
            sim = Simulation(**workflow_data)

            # Add to knowledge graph
            sim_uri = sim.to_graph(self.kg)
            workflow_uris.append(sim_uri)

            print(
                f"Added workflow {i+1}: connecting samples "
                f"{workflow_data.get('inputs', [])} to {workflow_data.get('outputs', [])}"
            )

        return workflow_uris

    def parse_activities(self, activity_data_list: List[Dict[str, Any]]) -> List[str]:
        """
        Parse activity data (transformations between samples).

        Parameters
        ----------
        activity_data_list : list of dict
            List of activity dictionaries

        Returns
        -------
        list of str
            List of activity URIs created

        Notes
        -----
        This method is not yet fully implemented.
        """
        # TODO: Implement activity parsing
        # This will be similar to workflow parsing but for Activity objects
        activity_uris = []
        print(
            f"Activity parsing not yet implemented. Found {len(activity_data_list)} activities."
        )
        return activity_uris

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
            - 'activity_uris' : list of created activity URIs

        Raises
        ------
        ValueError
            If file format is not supported (must be .yaml, .yml, or .json)
        TypeError
            If data type is not supported
        """
        # If data is a file path, read it first
        if isinstance(data, (str, Path)):
            filepath = Path(data)

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

        result = {"sample_map": {}, "workflow_uris": [], "activity_uris": []}

        # Parse samples first (they may be referenced by workflows/activities)
        if "computational_sample" in data:
            result["sample_map"] = self.parse_samples(data["computational_sample"])

        # Parse workflows
        if "workflow" in data:
            result["workflow_uris"] = self.parse_workflows(data["workflow"])

        # Parse activities
        if "activity" in data:
            result["activity_uris"] = self.parse_activities(data["activity"])

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


def parse_activity(
    activity_data: Dict[str, Any], graph: Optional[KnowledgeGraph] = None
) -> str:
    """
    Parse a single activity.

    Parameters
    ----------
    activity_data : dict
        Activity dictionary
    graph : KnowledgeGraph, optional
        Knowledge graph instance. If None, creates a new one.

    Returns
    -------
    str or None
        URI of the created activity
    """
    parser = WorkflowParser(kg=graph)
    uris = parser.parse_activities([activity_data])
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
