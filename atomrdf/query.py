import numpy as np
from atomrdf.build.bulk import _generate_atomic_sample_data
from atomrdf.datamodels.structure import AtomicScaleSample


def materials_project(
    api_key,
    chemical_system=None,
    material_ids=None,
    is_stable=True,
    conventional=True,
    graph=None,
):
    try:
        from mp_api.client import MPRester
    except ImportError:
        raise ImportError(
            "mp-api is not installed. Please install it for this functionality."
        )
    rest = {
        "use_document_model": False,
        "include_user_agent": True,
        "api_key": api_key,
    }
    if (chemical_system is None) and (material_ids is None):
        raise ValueError(
            "Please provide either a chemical system or a list of material ids"
        )

    with MPRester(**rest) as mpr:
        if chemical_system is not None:
            docs = mpr.materials.summary.search(
                chemsys=chemical_system, is_stable=is_stable
            )
        else:
            docs = mpr.materials.summary.search(material_ids=material_ids)

    # process docs
    structures = []

    for doc in docs:
        struct = doc["structure"]
        if conventional:
            aseatoms = struct.to_conventional().to_ase_atoms()
        else:
            aseatoms = struct.to_primitive().to_ase_atoms()

        symmetry = doc["symmetry"]

        if graph is not None:
            data = _generate_atomic_sample_data(aseatoms)
            sys = AtomicScaleSample(**data)
            sys.material.crystal_structure.spacegroup_symbol = symmetry["symbol"]
            sys.material.crystal_structure.spacegroup_number = symmetry["number"]
            sys.to_graph(graph)
            aseatoms.info["id"] = sys.id
            aseatoms.info["graph"] = graph

            # add energy
            # TODO
            # self.add_calculated_quantity(
            #    sys.sample, "EnergyPerAtom", doc["energy_per_atom"], unit="EV"
            # )
            structures.append(sys)
    if len(structures) == 1:
        return structures[0]
    else:
        return structures
