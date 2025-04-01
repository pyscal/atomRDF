from ase.io import write
from pyscal3.formats.ase import convert_snap

def _write_poscar_id(system, outfile):
    lines = []
    with open(outfile, "r") as fin:
        for line in fin:
            lines.append(line)
    lines[0] = system.sample.toPython() + "\n"
    with open(outfile, "w") as fout:
        for line in lines:
            fout.write(line)

def write(system, filename=None, 
            format='lammps-dump', 
            customkeys=None, customvals=None,
            compressed=False, timestep=0, species=None,  add_sample_id=True,
            copy_from=None, pseudo_files=None):
    """
    Write the structure to a file in the specified format.

    Parameters
    ----------
    outfile : str
        The path to the output file.
    format : str, optional
        The format of the output file. Defaults to 'lammps-dump'.
    customkeys : list, optional
        A list of custom keys to include in the output file. Defaults to None.
        Only valid if format is 'lammps-dump'.
    customvals : list, optional
        A list of custom values corresponding to the custom keys. Defaults to None.
        Only valid if format is 'lammps-dump'.
    compressed : bool, optional
        Whether to compress the output file. Defaults to False.
    timestep : int, optional
        The timestep value to include in the output file. Defaults to 0.
        Only valid if format is 'lammps-dump'.
    species : list, optional
        A list of species to include in the output file. Defaults to None.
        Only valid for ASE, if species is not specified.
    add_sample_id : bool, optional
        Whether to add a sample ID to the output file. Defaults to True.
        Only valid for poscar and quantum-espresso formats.
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
    if filename is None:
        outfile = f"structure.out"
    else:
        outfile = filename

    if format == "ase":
        asesys = convert_snap(system)
        if system.sample is not None:
            asesys.info["sample_id"] = system.sample
        return asesys
    
    elif format == "pyiron":
        from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
        asesys = convert_snap(system)
        pyironsys = ase_to_pyiron(asesys)
        if system.sample is not None:
            pyironsys.info["sample_id"] = system.sample            
        return pyironsys

    elif format == "poscar":
        asesys = convert_snap(system)
        write(outfile, asesys, format="vasp")
        if add_sample_id and (system.sample is not None):
            system.write_poscar_id(outfile)
    
    elif format == "lammps-dump":
        inputmethods.to_file(system, outfile, format='lammps-dump', customkeys=customkeys, customvals=customvals,
            compressed=compressed, timestep=timestep, species=species)
    
    elif format == "lammps-data":
        asesys = convert_snap(system)
        write(outfile, asesys, format='lammps-data', atom_style='atomic')
    
    elif format == "quantum-espresso":
        aio.write_espresso(system, filename, copy_from=copy_from, pseudo_files=pseudo_files)
    
    else:
        asesys = convert_snap(system)
        write(outfile, asesys, format=format)