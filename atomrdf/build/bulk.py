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
    get_metadata: bool = False,
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
    sdict["spacegroup_symbol"] = ap.get_spacegroup_symbol(atoms)
    sdict["spacegroup_number"] = ap.get_spacegroup_number(atoms)
    data = _generate_atomic_sample_data(atoms, sdict, repeat)

    if graph is not None:
        sample = AtomicScaleSample(**data)
        sample.to_graph(graph)
        atoms.info["id"] = sample.id
        atoms.info["graph"] = graph

    if get_metadata:
        return atoms, sdict
    return atoms


def _generate_atomic_sample_data(atoms, sdict=None, repeat=None):
    data = AtomicScaleSample.template()
    data["material"]["element_ratio"] = ap.get_chemical_composition(atoms)

    if sdict is not None:
        if "structure" in sdict.keys():
            data["material"]["crystal_structure"]["name"] = sdict["structure"]
        if "spacegroup_symbol" in sdict.keys():
            data["material"]["crystal_structure"]["spacegroup_symbol"] = sdict[
                "spacegroup_symbol"
            ]
        if "spacegroup_number" in sdict.keys():
            data["material"]["crystal_structure"]["spacegroup_number"] = sdict[
                "spacegroup_number"
            ]
        if "structure" in sdict.keys():
            data["material"]["crystal_structure"]["unit_cell"]["bravais_lattice"] = (
                ap.get_bravais_lattice(sdict["structure"])
            )
        if "a" in sdict.keys():
            if "b" not in sdict.keys():
                sdict["b"] = sdict["a"]
            if "c" not in sdict.keys():
                sdict["c"] = sdict["a"]
            data["material"]["crystal_structure"]["unit_cell"]["lattice_parameter"] = [
                sdict["a"],
                sdict["b"],
                sdict["c"],
            ]

    data["material"]["crystal_structure"]["unit_cell"][
        "angle"
    ] = atoms.get_cell_lengths_and_angles()[3:]

    data["simulation_cell"]["volume"]["value"] = ap.get_cell_volume(atoms)
    data["simulation_cell"]["number_of_atoms"] = ap.get_number_of_atoms(atoms)
    data["simulation_cell"]["length"] = ap.get_simulation_cell_length(atoms)
    data["simulation_cell"]["vector"] = ap.get_simulation_cell_vector(atoms)
    data["simulation_cell"]["angle"] = ap.get_simulation_cell_angle(atoms)

    if repeat is not None:
        if isinstance(repeat, int):
            data["simulation_cell"]["repetitions"] = (repeat, repeat, repeat)
        else:
            data["simulation_cell"]["repetitions"] = repeat

    data["atom_attribute"]["position"] = atoms.get_positions().tolist()
    data["atom_attribute"]["species"] = atoms.get_chemical_symbols()
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
    if isinstance(repeat, (int, list, tuple)):
        atoms = atoms.repeat(repeat)
    return atoms


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
            covera = ref.get("c/a") if xref == crystalstructure else np.sqrt(8 / 3)

    if covera is None and ref and (ref_c_a := ref.get("c/a")):
        covera = ref_c_a
        if c is None and a:
            sdict["c"] = covera * a

    sdict["b"] = sdict["a"] if sdict["b"] is None else sdict["b"]
    sdict["c"] = sdict["a"] if sdict["c"] is None else sdict["c"]

    return sdict
