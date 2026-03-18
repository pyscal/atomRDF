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
    final_structure = read("tmp.dump", format="lammps-dump-text")
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


def calculate_ev_curves(
    structure,
    pair_style,
    pair_coeff,
    vol_range=0.3,
    num_of_points=5,
    cores=1,
    e_tol=0,
    f_tol=0.0001,
    n_energy_steps=1e5,
    n_force_steps=1e6,
):
    """
    Calculate energy-volume curves for a structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input atomic structure
    pair_style : str
        LAMMPS pair style
    pair_coeff : str
        LAMMPS pair coefficients
    vol_range : float, optional
        Volume range for scaling (default: 0.3)
    num_of_points : int, optional
        Number of points in the curve (default: 5)
    cores : int, optional
        Number of cores for LAMMPS (default: 1)
    e_tol : float, optional
        Energy tolerance for minimization (default: 0)
    f_tol : float, optional
        Force tolerance for minimization (default: 0.0001)
    n_energy_steps : float, optional
        Maximum energy minimization steps (default: 1e5)
    n_force_steps : float, optional
        Maximum force minimization steps (default: 1e6)

    Returns
    -------
    dict
        Dictionary containing 'energy', 'volume', and 'bulk_modulus' arrays
    """
    from mendeleev import element
    import pandas as pd
    import numpy as np
    from scipy.optimize import curve_fit
    from ase.io import read, write
    from pylammpsmpi import LammpsLibrary

    volume_factors = np.linspace((1 - vol_range), (1.0 + vol_range), num_of_points)

    energies = []
    volumes = []
    initial_cell = structure.get_cell()
    initial_volume = structure.get_volume()

    for volume_factor in volume_factors:
        scaled_atoms = scale_atoms(structure, volume_factor)
        e, v = calculate_energy_rigid(
            scaled_atoms,
            pair_style,
            pair_coeff,
            cores=cores,
            e_tol=e_tol,
            f_tol=f_tol,
            n_energy_steps=n_energy_steps,
            n_force_steps=n_force_steps,
        )

        energies.append(e)
        volumes.append(v)
    V0, E0, B0, Bp = fit_bm(volumes, energies)
    v_fit = np.linspace(min(volumes), max(volumes), 100, endpoint=True)
    e_fit = birch_murnaghan_eval(v_fit, V0, E0, B0, Bp)
    bulk_modulus = B0 * 160.2176621
    datadict = {"energy": e_fit, "volume": v_fit, "bulk_modulus": bulk_modulus}
    return datadict


def run_md_nvt(
    structure,
    pair_style,
    pair_coeff,
    temperature=300.0,
    n_equilibration_steps=10000,
    n_production_steps=10000,
    timestep=0.002,
    cores=1,
    seed=12345,
):
    """
    Run an NVT (canonical ensemble) molecular dynamics simulation.

    The thermostat is a Nosé-Hoover chain (``fix nvt``).  The structure is
    first equilibrated for *n_equilibration_steps* steps, then statistics
    are collected over *n_production_steps* steps.  The time-averaged
    potential energy per atom and volume per atom are returned together with
    the final structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure.
    pair_style : str
        LAMMPS pair style string.
    pair_coeff : str
        LAMMPS pair coefficient string.
    temperature : float, optional
        Target temperature in Kelvin (default: 300.0).
    n_equilibration_steps : int, optional
        Number of MD steps for thermostat equilibration (default: 10000).
    n_production_steps : int, optional
        Number of MD steps over which averages are collected (default: 10000).
    timestep : float, optional
        MD timestep in picoseconds (default: 0.002).
    cores : int, optional
        Number of MPI cores for LAMMPS (default: 1).
    seed : int, optional
        Random seed for velocity initialisation (default: 12345).

    Returns
    -------
    tuple
        (final_structure, mean_energy_per_atom, mean_volume_per_atom)
    """
    from mendeleev import element
    import numpy as np
    from ase.io import read, write
    from pylammpsmpi import LammpsLibrary

    symbols, counts = np.unique(structure.get_chemical_symbols(), return_counts=True)
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
    pair_coeff = np.atleast_1d(pair_coeff)
    for pair in pair_coeff:
        lmp.command(f"pair_coeff {pair}")

    lmp.command(f"velocity all create {temperature} {seed} dist gaussian")
    lmp.command(f"fix ensemble all nvt temp {temperature} {temperature} $(100*dt)")
    lmp.command(f"timestep {timestep}")

    # Equilibration
    lmp.command(f"run {int(n_equilibration_steps)}")

    # Production: accumulate averages and write to file
    avg_file = "tmp_nvt_avg.txt"
    lmp.command("variable pe_atom equal pe/atoms")
    lmp.command("variable vol_atom equal vol/atoms")
    lmp.command(
        f"fix avg_pe all ave/time 1 {int(n_production_steps)} {int(n_production_steps)}"
        f" v_pe_atom v_vol_atom file {avg_file}"
    )
    lmp.command(f"run {int(n_production_steps)}")
    lmp.command("unfix avg_pe")
    lmp.command("unfix ensemble")

    dump_file = "tmp_nvt.dump"
    lmp.command(
        f"dump final all custom 1 {dump_file} id type mass x y z vx vy vz"
    )
    lmp.command("run 0")
    lmp.command("undump final")
    lmp.close()

    # Read averages from file — last data line has the final accumulated values
    with open(avg_file) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    vals = lines[-1].split()
    mean_energy_per_atom = float(vals[1])
    mean_volume_per_atom = float(vals[2])

    final_structure = read(dump_file, format="lammps-dump-text")
    return final_structure, mean_energy_per_atom, mean_volume_per_atom


def run_md_npt(
    structure,
    pair_style,
    pair_coeff,
    temperature=300.0,
    pressure=0.0,
    n_equilibration_steps=10000,
    n_production_steps=10000,
    timestep=0.002,
    cores=1,
    seed=12345,
):
    """
    Run an NPT (isothermal-isobaric ensemble) molecular dynamics simulation.

    Uses a Nosé-Hoover thermostat + barostat (``fix npt``) with isotropic
    pressure coupling.  Equilibration is followed by a production run from
    which time-averaged potential energy per atom and volume per atom are
    returned together with the final structure.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure.
    pair_style : str
        LAMMPS pair style string.
    pair_coeff : str
        LAMMPS pair coefficient string.
    temperature : float, optional
        Target temperature in Kelvin (default: 300.0).
    pressure : float, optional
        Target pressure in bar (default: 0.0).
    n_equilibration_steps : int, optional
        Number of MD steps for thermostat/barostat equilibration (default: 10000).
    n_production_steps : int, optional
        Number of MD steps over which averages are collected (default: 10000).
    timestep : float, optional
        MD timestep in picoseconds (default: 0.002).
    cores : int, optional
        Number of MPI cores for LAMMPS (default: 1).
    seed : int, optional
        Random seed for velocity initialisation (default: 12345).

    Returns
    -------
    tuple
        (final_structure, mean_energy_per_atom, mean_volume_per_atom)
    """
    from mendeleev import element
    import numpy as np
    from ase.io import read, write
    from pylammpsmpi import LammpsLibrary

    symbols, counts = np.unique(structure.get_chemical_symbols(), return_counts=True)
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
    pair_coeff = np.atleast_1d(pair_coeff)
    for pair in pair_coeff:
        lmp.command(f"pair_coeff {pair}")

    lmp.command(f"velocity all create {temperature} {seed} dist gaussian")
    lmp.command(
        f"fix ensemble all npt"
        f" temp {temperature} {temperature} $(100*dt)"
        f" iso {pressure} {pressure} $(1000*dt)"
    )
    lmp.command(f"timestep {timestep}")

    # Equilibration
    lmp.command(f"run {int(n_equilibration_steps)}")

    # Production: accumulate averages and write to file
    avg_file = "tmp_npt_avg.txt"
    lmp.command("variable pe_atom equal pe/atoms")
    lmp.command("variable vol_atom equal vol/atoms")
    lmp.command(
        f"fix avg_npt all ave/time 1 {int(n_production_steps)} {int(n_production_steps)}"
        f" v_pe_atom v_vol_atom file {avg_file}"
    )
    lmp.command(f"run {int(n_production_steps)}")
    lmp.command("unfix avg_npt")
    lmp.command("unfix ensemble")

    dump_file = "tmp_npt.dump"
    lmp.command(
        f"dump final all custom 1 {dump_file} id type mass x y z vx vy vz"
    )
    lmp.command("run 0")
    lmp.command("undump final")
    lmp.close()

    # Read averages from file — last data line has the final accumulated values
    with open(avg_file) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    vals = lines[-1].split()
    mean_energy_per_atom = float(vals[1])
    mean_volume_per_atom = float(vals[2])

    final_structure = read(dump_file, format="lammps-dump-text")
    return final_structure, mean_energy_per_atom, mean_volume_per_atom


def scale_atoms(structure, scale_factor):
    """
    Scale atomic structure uniformly.

    Parameters
    ----------
    structure : ase.Atoms
        Input structure
    scale_factor : float
        Scaling factor

    Returns
    -------
    ase.Atoms
        Scaled structure
    """
    scaled_atoms = structure.copy()
    scaled_atoms.set_cell(scaled_atoms.get_cell() * scale_factor, scale_atoms=True)
    return scaled_atoms


def fit_bm(vol, en):
    """
    Fit Birch-Murnaghan equation of state.

    Parameters
    ----------
    vol : array-like
        Volume values
    en : array-like
        Energy values

    Returns
    -------
    tuple
        (V0, E0, B0, Bp) - equilibrium volume, energy, bulk modulus, and derivative
    """
    import numpy as np
    from scipy.optimize import curve_fit

    a, b, c = np.polyfit(vol, en, 2)
    V0 = -b / (2 * a)
    E0 = a * V0**2 + b * V0 + c
    B0 = 2 * a * V0
    Bp = 4.0
    popt, pcov = curve_fit(birch_murnaghan_eval, vol, en, p0=[V0, E0, B0, Bp])
    return popt


def birch_murnaghan_eval(vol, V0, E0, B0, Bp):
    """
    Evaluate Birch-Murnaghan equation of state.

    Parameters
    ----------
    vol : array-like
        Volume values
    V0 : float
        Equilibrium volume
    E0 : float
        Equilibrium energy
    B0 : float
        Bulk modulus
    Bp : float
        Derivative of bulk modulus

    Returns
    -------
    array-like
        Energy values
    """
    eta = (vol / V0) ** (1.0 / 3.0)
    E = E0 + 9.0 * B0 * V0 / 16.0 * (eta**2 - 1.0) ** 2 * (
        6.0 + Bp * (eta**2 - 1.0) - 4.0 * eta**2
    )
    return E
