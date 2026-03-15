def calculate_energy_relax(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
):
    """
    Relax structure with box relaxation.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    cores : int, optional
        Number of cores (default: 1)
    e_tol : float, optional
        Energy tolerance (default: 0)
    f_tol : float, optional
        Force tolerance (default: 0.0001)
    n_energy_steps : float, optional
        Max energy steps (default: 1e5)
    n_force_steps : float, optional
        Max force steps (default: 1e6)

    Returns
    -------
    tuple
        (final_structure, ecoh, vol) - relaxed structure, cohesive energy, volume
    """
    from mendeleev import element
    import pandas as pd
    import numpy as np
    from scipy.optimize import curve_fit
    from ase.io import read, write
    from pylammpsmpi import LammpsLibrary

    symbols, counts = list(
        np.unique(structure.get_chemical_symbols(), return_counts=True)
    )
    masses = [element(symbol).mass for symbol in symbols]
    write("tmp.data", structure, format="lammps-data")

    initial_cell = structure.get_cell()
    initial_volume = structure.get_volume()

    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    for x in range(len(masses)):
        lmp.command(f"mass {x+1} {masses[x]}")

    lmp.command(f"pair_style {pair_style}")

    pair_coeff = np.atleast_1d(pair_coeff)
    for pair in pair_coeff:
        lmp.command(f"pair_coeff {pair}")

    lmp.command("fix ensemble all box/relax aniso 0.0")
    lmp.command(f"minimize {e_tol} {f_tol} {int(n_energy_steps)} {int(n_force_steps)}")

    lmp.command("run 0")
    ecoh = lmp.pe / lmp.natoms
    final_volume = lmp.vol
    vol = lmp.vol / lmp.natoms

    filename = "tmp.dump"
    lmp.command(
        "dump              2x all custom 1 %s id type mass x y z vx vy vz" % (filename)
    )
    lmp.command("run               0")
    lmp.command("undump            2x")

    lmp.close()
    return final_structure, ecoh, vol


def calculate_energy_rigid(
    structure,
    pair_style,
    pair_coeff,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
):
    """
    Calculate energy of a structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    cores : int, optional
        Number of cores (default: 1)
    e_tol : float, optional
        Energy tolerance (default: 0)
    f_tol : float, optional
        Force tolerance (default: 0.0001)
    n_energy_steps : float, optional
        Max energy steps (default: 1e5)
    n_force_steps : float, optional
        Max force steps (default: 1e6)

    Returns
    -------
    tuple
        (ecoh, vol) - cohesive energy and volume per atom
    """
    from mendeleev import element
    import pandas as pd
    import numpy as np
    from scipy.optimize import curve_fit
    from ase.io import read, write
    from pylammpsmpi import LammpsLibrary

    symbols, counts = list(
        np.unique(structure.get_chemical_symbols(), return_counts=True)
    )
    masses = [element(symbol).mass for symbol in symbols]

    write("tmp.data", structure, format="lammps-data")
    lmp = LammpsLibrary(cores=cores)
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("read_data tmp.data")
    for x in range(len(masses)):
        lmp.command(f"mass {x+1} {masses[x]}")

    lmp.command(f"pair_style {pair_style}")
    lmp.command(f"pair_coeff {pair_coeff}")
    lmp.command("run 0")
    ecoh = lmp.pe / lmp.natoms
    vol = lmp.vol / lmp.natoms
    return ecoh, vol
