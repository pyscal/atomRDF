"""
Integration test for workflow reconstruction.

Builds a small KG from a YAML, then calls reconstruct_workflow to
verify that structure files, input files, and metadata are written.
"""

import os
import shutil
import tempfile
import yaml

from atomrdf import KnowledgeGraph
from atomrdf.io.workflow_parser import WorkflowParser


# -- A minimal workflow YAML for testing --
TEST_YAML = {
    "computational_sample": [
        {
            "id": "pristine_Al",
            "simulation_cell": {
                "number_of_atoms": 4,
                "length": [4.05, 4.05, 4.05],
                "angle": [90, 90, 90],
                "vector": [[4.05, 0, 0], [0, 4.05, 0], [0, 0, 4.05]],
            },
            "material": {
                "crystal_structure": {
                    "name": "FCC",
                    "spacegroup_symbol": "Fm-3m",
                    "spacegroup_number": 225,
                    "unit_cell": {
                        "bravais_lattice": "cF",
                        "lattice_parameter": [4.05, 4.05, 4.05],
                        "angle": [90, 90, 90],
                    },
                },
                "element_ratio": {"Al": 1.0},
            },
            "atom_attribute": {
                "species": [
                    {
                        "chemical_symbol": "Al",
                        "positions": [
                            [0, 0, 0],
                            [0.5, 0.5, 0],
                            [0.5, 0, 0.5],
                            [0, 0.5, 0.5],
                        ],
                    }
                ],
            },
        },
        {
            "id": "relaxed_Al",
            "simulation_cell": {
                "number_of_atoms": 4,
                "length": [4.04, 4.04, 4.04],
                "angle": [90, 90, 90],
                "vector": [[4.04, 0, 0], [0, 4.04, 0], [0, 0, 4.04]],
            },
            "material": {
                "crystal_structure": {
                    "name": "FCC",
                    "spacegroup_symbol": "Fm-3m",
                    "spacegroup_number": 225,
                    "unit_cell": {
                        "bravais_lattice": "cF",
                        "lattice_parameter": [4.04, 4.04, 4.04],
                        "angle": [90, 90, 90],
                    },
                },
                "element_ratio": {"Al": 1.0},
            },
            "atom_attribute": {
                "species": [
                    {
                        "chemical_symbol": "Al",
                        "positions": [
                            [0, 0, 0],
                            [0.5, 0.5, 0],
                            [0.5, 0, 0.5],
                            [0, 0.5, 0.5],
                        ],
                    }
                ],
            },
        },
    ],
    "workflow": [
        {
            "method": "DensityFunctionalTheory",
            "software": [{"label": "VASP", "uri": "https://www.vasp.at"}],
            "degrees_of_freedom": ["AtomicPositionRelaxation", "CellVolumeRelaxation"],
            "xc_functional": "GGA",
            "input_sample": ["pristine_Al"],
            "output_sample": ["relaxed_Al"],
            "input_parameter": [
                {"label": "ENCUT", "value": 400},
                {"label": "ISMEAR", "value": 1},
                {"label": "SIGMA", "value": 0.2},
            ],
            "calculated_property": [
                {
                    "label": "TotalEnergy",
                    "value": -14.123,
                    "unit": "EV",
                    "associate_to_sample": ["relaxed_Al"],
                },
            ],
        },
    ],
}


def test_single_step_vasp():
    """Reconstruct a single-step DFT workflow."""
    tmpdir = tempfile.mkdtemp(prefix="atomrdf_test_")
    try:
        struct_store = os.path.join(tmpdir, "structures")
        kg = KnowledgeGraph(structure_store=struct_store)
        parser = WorkflowParser(kg=kg)

        yaml_path = os.path.join(tmpdir, "test.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump(TEST_YAML, f)

        result = parser.parse(yaml_path)
        assert len(result["workflow_uris"]) == 1, "Expected 1 workflow"

        # Now reconstruct
        out_dir = os.path.join(tmpdir, "reconstructed")
        wf_uri = result["workflow_uris"][0]
        kg.reconstruct_workflow(wf_uri, out_dir)

        # Check outputs
        assert os.path.isdir(out_dir)
        assert os.path.isfile(os.path.join(out_dir, "INCAR")), "Missing INCAR"
        assert os.path.isfile(os.path.join(out_dir, "README.md")), "Missing README"
        assert os.path.isfile(
            os.path.join(out_dir, "software_requirements.txt")
        ), "Missing software_requirements.txt"

        # Check INCAR contents
        with open(os.path.join(out_dir, "INCAR")) as f:
            incar = f.read()
        assert "ENCUT" in incar
        assert "ISMEAR" in incar
        assert "SIGMA" in incar

        # Check README contents
        with open(os.path.join(out_dir, "README.md")) as f:
            readme = f.read()
        assert "DensityFunctionalTheory" in readme
        assert "VASP" in readme

        # Check software_requirements
        with open(os.path.join(out_dir, "software_requirements.txt")) as f:
            reqs = f.read()
        assert "VASP" in reqs

        # Check at least one structure file exists
        structure_files = [
            f
            for f in os.listdir(out_dir)
            if f.startswith("POSCAR") or f.endswith(".lmp")
        ]
        assert len(structure_files) > 0, "No structure files written"

        print("PASS: test_single_step_vasp")
    finally:
        shutil.rmtree(tmpdir)


def test_reconstruct_by_sample():
    """Reconstruct via sample ID."""
    tmpdir = tempfile.mkdtemp(prefix="atomrdf_test_")
    try:
        struct_store = os.path.join(tmpdir, "structures")
        kg = KnowledgeGraph(structure_store=struct_store)
        parser = WorkflowParser(kg=kg)

        yaml_path = os.path.join(tmpdir, "test.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump(TEST_YAML, f)

        result = parser.parse(yaml_path)

        # Get the output sample URI
        output_sample = result["sample_map"]["relaxed_Al"]

        out_dir = os.path.join(tmpdir, "reconstructed_by_sample")
        kg.reconstruct_workflow_by_sample(output_sample, out_dir)

        assert os.path.isdir(out_dir)
        assert os.path.isfile(os.path.join(out_dir, "INCAR"))
        assert os.path.isfile(os.path.join(out_dir, "README.md"))

        print("PASS: test_reconstruct_by_sample")
    finally:
        shutil.rmtree(tmpdir)


# -- Multi-step YAML: pristine → defect (operation) → relax (workflow 1) → static (workflow 2) --
MULTI_STEP_YAML = {
    "computational_sample": [
        {
            "id": "pristine_Cu",
            "simulation_cell": {
                "number_of_atoms": 32,
                "length": [7.22, 7.22, 7.22],
                "angle": [90, 90, 90],
                "vector": [[7.22, 0, 0], [0, 7.22, 0], [0, 0, 7.22]],
            },
            "material": {
                "crystal_structure": {
                    "name": "FCC",
                    "spacegroup_symbol": "Fm-3m",
                    "spacegroup_number": 225,
                    "unit_cell": {
                        "bravais_lattice": "cF",
                        "lattice_parameter": [3.61, 3.61, 3.61],
                        "angle": [90, 90, 90],
                    },
                },
                "element_ratio": {"Cu": 1.0},
            },
            "atom_attribute": {
                "species": [{"chemical_symbol": "Cu", "positions": [[0, 0, 0]]}],
            },
        },
        {
            "id": "vacancy_Cu_unrelaxed",
            "simulation_cell": {
                "number_of_atoms": 31,
                "length": [7.22, 7.22, 7.22],
                "angle": [90, 90, 90],
                "vector": [[7.22, 0, 0], [0, 7.22, 0], [0, 0, 7.22]],
            },
            "material": {
                "element_ratio": {"Cu": 1.0},
            },
            "atom_attribute": {
                "species": [{"chemical_symbol": "Cu", "positions": [[0, 0, 0]]}],
            },
        },
        {
            "id": "vacancy_Cu_relaxed",
            "simulation_cell": {
                "number_of_atoms": 31,
                "length": [7.20, 7.20, 7.20],
                "angle": [90, 90, 90],
                "vector": [[7.20, 0, 0], [0, 7.20, 0], [0, 0, 7.20]],
            },
            "material": {
                "element_ratio": {"Cu": 1.0},
            },
            "atom_attribute": {
                "species": [{"chemical_symbol": "Cu", "positions": [[0, 0, 0]]}],
            },
        },
        {
            "id": "vacancy_Cu_static",
            "simulation_cell": {
                "number_of_atoms": 31,
                "length": [7.20, 7.20, 7.20],
                "angle": [90, 90, 90],
                "vector": [[7.20, 0, 0], [0, 7.20, 0], [0, 0, 7.20]],
            },
            "material": {
                "element_ratio": {"Cu": 1.0},
            },
            "atom_attribute": {
                "species": [{"chemical_symbol": "Cu", "positions": [[0, 0, 0]]}],
            },
        },
    ],
    "operation": [
        {
            "method": "DeleteAtom",
            "input_sample": "pristine_Cu",
            "output_sample": "vacancy_Cu_unrelaxed",
        },
    ],
    "workflow": [
        {
            "method": "DensityFunctionalTheory",
            "software": [{"label": "VASP", "uri": "https://www.vasp.at"}],
            "degrees_of_freedom": ["AtomicPositionRelaxation", "CellVolumeRelaxation"],
            "xc_functional": "GGA",
            "input_sample": ["vacancy_Cu_unrelaxed"],
            "output_sample": ["vacancy_Cu_relaxed"],
            "input_parameter": [
                {"label": "ENCUT", "value": 520},
                {"label": "EDIFF", "value": 1e-6},
            ],
            "calculated_property": [
                {
                    "label": "TotalEnergy",
                    "value": -110.5,
                    "unit": "EV",
                    "associate_to_sample": ["vacancy_Cu_relaxed"],
                },
            ],
        },
        {
            "method": "DensityFunctionalTheory",
            "software": [{"label": "VASP", "uri": "https://www.vasp.at"}],
            "degrees_of_freedom": [],
            "xc_functional": "GGA",
            "input_sample": ["vacancy_Cu_relaxed"],
            "output_sample": ["vacancy_Cu_static"],
            "input_parameter": [
                {"label": "ENCUT", "value": 520},
                {"label": "EDIFF", "value": 1e-8},
                {"label": "NSW", "value": 0},
            ],
            "calculated_property": [
                {
                    "label": "TotalEnergy",
                    "value": -110.4,
                    "unit": "EV",
                    "associate_to_sample": ["vacancy_Cu_static"],
                },
            ],
        },
    ],
}


def test_multistep_workflow():
    """Reconstruct a multi-step workflow (relax → static)."""
    tmpdir = tempfile.mkdtemp(prefix="atomrdf_test_")
    try:
        struct_store = os.path.join(tmpdir, "structures")
        kg = KnowledgeGraph(structure_store=struct_store)
        parser = WorkflowParser(kg=kg)

        yaml_path = os.path.join(tmpdir, "test_multi.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump(MULTI_STEP_YAML, f)

        result = parser.parse(yaml_path)
        assert len(result["workflow_uris"]) == 2

        # Reconstruct using the LAST workflow (static)
        wf_uri = result["workflow_uris"][1]
        out_dir = os.path.join(tmpdir, "reconstructed_multi")
        kg.reconstruct_workflow(wf_uri, out_dir)

        assert os.path.isdir(out_dir)

        # Should have step directories
        has_steps = os.path.isdir(os.path.join(out_dir, "step_1"))
        has_readme = os.path.isfile(os.path.join(out_dir, "README.md"))
        has_postprocess = os.path.isfile(os.path.join(out_dir, "postprocess.py"))
        has_reqs = os.path.isfile(os.path.join(out_dir, "software_requirements.txt"))

        assert has_readme, "Missing top-level README.md"
        assert has_reqs, "Missing software_requirements.txt"

        # Check software requirements
        with open(os.path.join(out_dir, "software_requirements.txt")) as f:
            reqs = f.read()
        assert "VASP" in reqs

        print("PASS: test_multistep_workflow")
    finally:
        shutil.rmtree(tmpdir)


def test_lammps_workflow():
    """Reconstruct a LAMMPS workflow."""
    lammps_yaml = {
        "computational_sample": [
            {
                "id": "input_Fe",
                "simulation_cell": {
                    "number_of_atoms": 2,
                    "length": [2.87, 2.87, 2.87],
                    "angle": [90, 90, 90],
                    "vector": [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]],
                },
                "material": {
                    "crystal_structure": {
                        "name": "BCC",
                        "spacegroup_symbol": "Im-3m",
                        "spacegroup_number": 229,
                        "unit_cell": {
                            "bravais_lattice": "cI",
                            "lattice_parameter": [2.87, 2.87, 2.87],
                            "angle": [90, 90, 90],
                        },
                    },
                    "element_ratio": {"Fe": 1.0},
                },
                "atom_attribute": {
                    "species": [
                        {
                            "chemical_symbol": "Fe",
                            "positions": [[0, 0, 0], [0.5, 0.5, 0.5]],
                        }
                    ],
                },
            },
            {
                "id": "relaxed_Fe",
                "simulation_cell": {
                    "number_of_atoms": 2,
                    "length": [2.86, 2.86, 2.86],
                    "angle": [90, 90, 90],
                    "vector": [[2.86, 0, 0], [0, 2.86, 0], [0, 0, 2.86]],
                },
                "material": {
                    "element_ratio": {"Fe": 1.0},
                },
                "atom_attribute": {
                    "species": [
                        {
                            "chemical_symbol": "Fe",
                            "positions": [[0, 0, 0], [0.5, 0.5, 0.5]],
                        }
                    ],
                },
            },
        ],
        "workflow": [
            {
                "method": "MolecularStatics",
                "software": [{"label": "LAMMPS", "uri": "https://www.lammps.org"}],
                "degrees_of_freedom": ["AtomicPositionRelaxation"],
                "input_sample": ["input_Fe"],
                "output_sample": ["relaxed_Fe"],
                "interatomic_potential": {
                    "uri": "https://doi.org/10.xxx",
                    "potential_type": "EmbeddedAtomModel",
                    "label": "Fe_EAM_2020",
                },
                "input_parameter": [
                    {"label": "minimize_etol", "value": 1e-10},
                ],
                "calculated_property": [
                    {
                        "label": "TotalEnergy",
                        "value": -8.5,
                        "unit": "EV",
                        "associate_to_sample": ["relaxed_Fe"],
                    },
                ],
            },
        ],
    }
    tmpdir = tempfile.mkdtemp(prefix="atomrdf_test_")
    try:
        struct_store = os.path.join(tmpdir, "structures")
        kg = KnowledgeGraph(structure_store=struct_store)
        parser = WorkflowParser(kg=kg)

        yaml_path = os.path.join(tmpdir, "lammps.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump(lammps_yaml, f)

        result = parser.parse(yaml_path)
        wf_uri = result["workflow_uris"][0]

        out_dir = os.path.join(tmpdir, "reconstructed_lammps")
        kg.reconstruct_workflow(wf_uri, out_dir)

        assert os.path.isdir(out_dir)
        assert os.path.isfile(os.path.join(out_dir, "in.lammps")), "Missing in.lammps"
        assert os.path.isfile(os.path.join(out_dir, "README.md")), "Missing README"

        # Check LAMMPS input
        with open(os.path.join(out_dir, "in.lammps")) as f:
            inp = f.read()
        assert "read_data" in inp
        assert "minimize" in inp

        # Check README for LAMMPS details
        with open(os.path.join(out_dir, "README.md")) as f:
            readme = f.read()
        assert "MolecularStatics" in readme
        assert "LAMMPS" in readme

        print("PASS: test_lammps_workflow")
    finally:
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_single_step_vasp()
    test_reconstruct_by_sample()
    test_multistep_workflow()
    test_lammps_workflow()
    print("\nAll tests passed!")
