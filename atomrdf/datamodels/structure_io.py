from ase import Atoms
from ase.io import write as ase_write
from ase.io.espresso import read_fortran_namelist
import numpy as np
import os
import warnings
import mendeleev


def sample_to_ase(sample):
    """
    Convert a CMSO AtomicScaleSample to an ASE Atoms object
    """
    # first get the cell
    cell = sample.simulation_cell.vector

    # now read atom attributes - positions and species
    positions = sample.atom_attribute.position
    species = sample.atom_attribute.species

    # create the ASE Atoms object
    ase_atoms = Atoms(positions=positions, symbols=species, cell=cell)
    if sample.id:
        ase_atoms.info["id"] = sample.id
    return ase_atoms


def _convert_tab_to_dict(tab):
    keywords = [
        "ATOMIC_SPECIES",
        "ATOMIC_POSITIONS",
        "K_POINTS",
        "CELL_PARAMETERS",
        "OCCUPATIONS",
        "CONSTRAINTS",
        "ATOMIC_VELOCITIES",
        "ATOMIC_FORCES",
        "ADDITIONAL_K_POINTS",
        "SOLVENTS",
        "HUBBARD",
    ]

    tabdict = {}
    for line in tab:
        firstword = line.split()[0]
        secondword = " ".join(line.split()[1:])

        if firstword in keywords:
            tabdict[firstword] = {}
            tabdict[firstword]["value"] = []
            tabdict[firstword]["extra"] = secondword
            tabarr = tabdict[firstword]["value"]
        else:
            tabarr.append(line.strip())
    return tabdict


def _write_espresso(atoms, inputfile, copy_from=None, pseudo_files=None):
    data = None
    tab = None

    if copy_from is not None:
        if os.path.exists(copy_from):
            try:
                with open(copy_from, "r") as fin:
                    data, tab = read_fortran_namelist(fin)
            except:
                warnings.warn(
                    f"Error reading {copy_from}, a clean file will be written"
                )
                copy = True

    if tab is not None:
        tab = _convert_tab_to_dict(tab)
    else:
        tab = {}

    if data is None:
        data = {}
        data["system"] = {}
        data["control"] = {}

    tab["CELL_PARAMETERS"] = {}
    tab["CELL_PARAMETERS"]["extra"] = "angstrom"
    tab["CELL_PARAMETERS"]["value"] = []

    for vec in s.box:
        tab["CELL_PARAMETERS"]["value"].append(" ".join([str(x) for x in vec]))

    cds = atoms.get_scaled_positions()
    species = atoms.get_chemical_symbols()

    unique_species = np.unique(species)
    if pseudo_files is not None:
        if not len(pseudo_files) == len(unique_species):
            raise ValueError(
                "Number of pseudo files must match number of unique species"
            )
        pseudo_dirs = [
            os.path.dirname(os.path.abspath(pseudo_file))
            for pseudo_file in pseudo_files
        ]
        if not len(np.unique(pseudo_dirs)) == 1:
            raise ValueError("All pseudo files must be in the same directory")
        data["control"]["pseudo_dir"] = pseudo_dirs[0]
    else:
        pseudo_files = ["None" for x in range(len(unique_species))]

    tab["ATOMIC_SPECIES"] = {}
    tab["ATOMIC_SPECIES"]["extra"] = ""
    tab["ATOMIC_SPECIES"]["value"] = []

    for count, us in enumerate(unique_species):
        chem = mendeleev.element(us)
        tab["ATOMIC_SPECIES"]["value"].append(
            f"{us} {chem.atomic_weight} {os.path.basename(pseudo_files[count])}"
        )

    tab["ATOMIC_POSITIONS"] = {}
    tab["ATOMIC_POSITIONS"]["extra"] = "crystal"
    tab["ATOMIC_POSITIONS"]["value"] = []

    for cd, sp in zip(cds, species):
        tab["ATOMIC_POSITIONS"]["value"].append(f"{sp} {cd[0]} {cd[1]} {cd[2]}")

    data["system"]["ibrav"] = 0
    data["system"]["nat"] = len(species)
    data["system"]["ntyp"] = len(unique_species)

    with open(inputfile, "w") as fout:
        if "id" in atoms.info:
            fout.write(f"! {atoms.info['id']}\n\n")

        for key, val in data.items():
            fout.write(f"&{key.upper()}\n")
            for k, v in val.items():
                if isinstance(v, str):
                    fout.write(f"   {k} = '{v}',\n")
                else:
                    fout.write(f"   {k} = {v},\n")
            fout.write("/\n")
            fout.write("\n")

        for key, val in tab.items():
            fout.write(f'{key} {val["extra"]}\n')
            fout.write("\n")
            for v in val["value"]:
                fout.write(v)
                fout.write("\n")
            fout.write("\n")


def _write_poscar_id(atoms, outfile):
    if "id" in atoms.info:
        lines = []
        with open(outfile, "r") as fin:
            for line in fin:
                lines.append(line)
        lines[0] = atoms.info["id"] + "\n"
        with open(outfile, "w") as fout:
            for line in lines:
                fout.write(line)


def write(
    system,
    outfile,
    format,
    copy_from=None,
    pseudo_files=None,
):
    """
    Write the structure to a file in the specified format.

    Parameters
    ----------
    outfile : str
        The path to the output file.
    format : str, optional
        The format of the output file. Defaults to 'lammps-dump'.
    copy_from : str, optional
        If provided, input options for quantum-espresso format will be copied from
        the given file. Structure specific information will be replaced.
        Note that the validity of input file is not checked.
    pseudo_files : list, optional
        if provided, add the pseudopotential filenames to file.
        Should be in alphabetical order of chemical species symbols.

    Returns
    -------
    None
    """
    if format == "poscar":
        ase_write(outfile, system, format="vasp")
        _write_poscar_id(system, outfile)

    elif format == "lammps-dump":
        ase_write(outfile, system, format="lammps-dump-text")

    elif format == "lammps-data":
        ase_write(outfile, system, format="lammps-data", atom_style="atomic")

    elif format == "quantum-espresso":
        _write_espresso(system, outfile, copy_from=copy_from, pseudo_files=pseudo_files)

    else:
        ase_write(outfile, system, format=format)
