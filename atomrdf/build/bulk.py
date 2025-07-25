import numpy as np

from ase.data import atomic_numbers, chemical_symbols, reference_states
from ase.atoms import Atoms
from ase.build import bulk as asebulk

from atomrdf.datamodels.structure import AtomicScaleSample
import atomrdf.properties as ap


def bulk(
    name: str,
    crystalstructure: str = None,
    a: float = None,
    b: float = None,
    c: float = None,
    *,
    alpha: float = None,
    covera: float = None,
    u: float = None,
    orthorhombic: bool = False,
    cubic: bool = False,
    basis=None,
    repeat: int = 1,
    graph=None,
) -> tuple[Atoms, AtomicScaleSample]:

    sdict = _compute_structure_metadata(name, crystalstructure, a, b, c, covera)
    atoms = _create_atoms(
        name,
        crystalstructure,
        sdict["a"],
        sdict["b"],
        sdict["c"],
        alpha,
        covera,
        u,
        orthorhombic,
        cubic,
        basis,
        repeat,
    )
    data = _generate_atomic_sample_data(atoms, sdict, repeat)
    if graph is not None:
        sample = AtomicScaleSample(**data)
        sample.to_graph(graph)
        atoms.info["id"] = sample.id
        atoms.info["graph"] = graph
    return atoms


def _generate_atomic_sample_data(atoms, sdict, repeat):
    data = AtomicScaleSample.template()
    data["material"]["element_ratio"]["value"] = ap.get_chemical_composition(atoms)
    data["material"]["crystal_structure"]["name"]["value"] = sdict["structure"]
    data["material"]["crystal_structure"]["spacegroup_symbol"]["value"] = (
        ap.get_spacegroup_symbol(atoms)
    )
    data["material"]["crystal_structure"]["spacegroup_number"]["value"] = (
        ap.get_spacegroup_number(atoms)
    )
    data["material"]["crystal_structure"]["unit_cell"]["bravais_lattice"]["value"] = (
        ap.get_bravais_lattice(sdict["structure"])
    )
    data["material"]["crystal_structure"]["unit_cell"]["lattice_parameter"]["value"] = [
        sdict["a"],
        sdict["b"],
        sdict["c"],
    ]
    data["material"]["crystal_structure"]["unit_cell"]["angle"][
        "value"
    ] = atoms.get_cell_lengths_and_angles()[3:]

    data["simulation_cell"]["volume"]["value"] = ap.get_cell_volume(atoms)
    data["simulation_cell"]["number_of_atoms"]["value"] = ap.get_number_of_atoms(atoms)
    data["simulation_cell"]["length"]["value"] = ap.get_simulation_cell_length(atoms)
    data["simulation_cell"]["vector"]["value"] = ap.get_simulation_cell_vector(atoms)
    data["simulation_cell"]["angle"]["value"] = ap.get_simulation_cell_angle(atoms)
    if isinstance(repeat, int):
        data["simulation_cell"]["repetitions"]["value"] = (repeat, repeat, repeat)
    else:
        data["simulation_cell"]["repetitions"]["value"] = repeat

    data["atom_attribute"]["position"]["value"] = atoms.get_positions().tolist()
    data["atom_attribute"]["species"]["value"] = atoms.get_chemical_symbols()
    return data


def _create_atoms(
    name,
    crystalstructure,
    a,
    b,
    c,
    alpha,
    covera,
    u,
    orthorhombic,
    cubic,
    basis,
    repeat,
):
    atoms = asebulk(
        name,
        crystalstructure=crystalstructure,
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        covera=covera,
        u=u,
        orthorhombic=orthorhombic,
        cubic=cubic,
        basis=basis,
    )
    return atoms.repeat((repeat,) * 3 if isinstance(repeat, int) else repeat)


def _compute_structure_metadata(name, crystalstructure, a, b, c, covera):
    sdict = {"a": a, "b": b, "c": c, "covera": covera}
    atomic_number = atomic_numbers.get(name)
    ref = reference_states[atomic_number]

    xref = None
    if ref:
        xref = ref.get("symmetry")
        if xref and name in chemical_symbols:
            sdict["structure"] = xref

    if crystalstructure:
        sdict["structure"] = crystalstructure

    if a is None and ref and "a" in ref:
        sdict["a"] = ref["a"]

    if b is None and ref and (bovera := ref.get("b/a")) and a:
        sdict["b"] = bovera * a

    if crystalstructure in ["hcp", "wurtzite"]:
        if c:
            covera = c / a
        elif covera is None:
            covera = ref.get("c/a") if xref == crystalstructure else sqrt(8 / 3)

    if covera is None and ref and (ref_c_a := ref.get("c/a")):
        covera = ref_c_a
        if c is None and a:
            sdict["c"] = covera * a

    sdict["b"] = sdict["a"] if sdict["b"] is None else sdict["b"]
    sdict["c"] = sdict["a"] if sdict["c"] is None else sdict["c"]

    return sdict
