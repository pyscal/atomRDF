from ase.io import write as ase_write
from pyscal3.formats.ase import convert_snap
import mendeleev
import numpy as np
from ase.io.espresso import read_fortran_namelist
import os
import warnings

def _convert_tab_to_dict(tab):
    keywords = ["ATOMIC_SPECIES",
            "ATOMIC_POSITIONS",
            "K_POINTS",
            "CELL_PARAMETERS",
            "OCCUPATIONS",
            "CONSTRAINTS",
            "ATOMIC_VELOCITIES",
            "ATOMIC_FORCES",
            "ADDITIONAL_K_POINTS",
            "SOLVENTS",
            "HUBBARD"]

    tabdict = {}
    for line in tab:
        firstword = line.split()[0]
        secondword = " ".join(line.split()[1:]) 

        if firstword in keywords:
            tabdict[firstword] = {}
            tabdict[firstword]['value'] = []
            tabdict[firstword]['extra'] = secondword
            tabarr = tabdict[firstword]['value']
        else:
            tabarr.append(line.strip())
    return tabdict

def _write_espresso(s, inputfile, copy_from=None, pseudo_files=None):
    data = None
    tab = None

    if copy_from is not None:
        if os.path.exists(copy_from):
            try:
                with open(copy_from, 'r') as fin:
                    data, tab = read_fortran_namelist(fin)
            except:
                warnings.warn(f'Error reading {copy_from}, a clean file will be written')
                copy=True
    
    if tab is not None:
        tab = _convert_tab_to_dict(tab)
    else:
        tab = {}
    
    if data is None:
        data = {}
        data['system'] = {}
        data['control'] = {}
    
    tab['CELL_PARAMETERS'] = {}
    tab['CELL_PARAMETERS']['extra'] = 'angstrom'
    tab['CELL_PARAMETERS']['value'] = []

    for vec in s.box:
        tab['CELL_PARAMETERS']['value'].append(' '.join([str(x) for x in vec])) 

    cds = s.direct_coordinates
    species = s.atoms.species

    unique_species = np.unique(species)
    if pseudo_files is not None:
        if not len(pseudo_files) == len(unique_species):
            raise ValueError('Number of pseudo files must match number of unique species')
        pseudo_dirs = [os.path.dirname(os.path.abspath(pseudo_file)) for pseudo_file in pseudo_files]
        if not len(np.unique(pseudo_dirs)) == 1:
            raise ValueError('All pseudo files must be in the same directory')
        data['control']['pseudo_dir'] = pseudo_dirs[0]
    else:
        pseudo_files = ['None' for x in range(len(unique_species))]

    tab['ATOMIC_SPECIES'] = {}
    tab['ATOMIC_SPECIES']['extra'] = ''
    tab['ATOMIC_SPECIES']['value'] = []

    for count, us in enumerate(unique_species):
        chem = mendeleev.element(us)
        tab['ATOMIC_SPECIES']['value'].append(f'{us} {chem.atomic_weight} {os.path.basename(pseudo_files[count])}')

    tab['ATOMIC_POSITIONS'] = {}
    tab['ATOMIC_POSITIONS']['extra'] = 'crystal'
    tab['ATOMIC_POSITIONS']['value'] = []

    for cd, sp in zip(cds, species):
        tab['ATOMIC_POSITIONS']['value'].append(f'{sp} {cd[0]} {cd[1]} {cd[2]}')

    data['system']['ibrav'] = 0
    data['system']['nat'] = len(species)
    data['system']['ntyp'] = len(unique_species)

    with open(inputfile, 'w') as fout:
        if s.sample is not None:
            fout.write(f'! {s.sample.toPython()}\n\n')
        
        for key, val in data.items():
            fout.write(f'&{key.upper()}\n')
            for k, v in val.items():
                if isinstance(v, str):
                    fout.write(f'   {k} = \'{v}\',\n')
                else:
                    fout.write(f'   {k} = {v},\n')
            fout.write('/\n')
            fout.write('\n')

        for key, val in tab.items():
            fout.write(f'{key} {val["extra"]}\n')
            fout.write('\n')
            for v in val['value']:
                fout.write(v)
                fout.write('\n')
            fout.write('\n')   

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
        ase_write(outfile, asesys, format="vasp")
        if add_sample_id and (system.sample is not None):
            _write_poscar_id(system, outfile)
    
    elif format == "lammps-dump":
        inputmethods.to_file(system, outfile, format='lammps-dump', customkeys=customkeys, customvals=customvals,
            compressed=compressed, timestep=timestep, species=species)
    
    elif format == "lammps-data":
        asesys = convert_snap(system)
        ase_write(outfile, asesys, format='lammps-data', atom_style='atomic')
    
    elif format == "quantum-espresso":
        _write_espresso(system, filename, copy_from=copy_from, pseudo_files=pseudo_files)
    
    else:
        asesys = convert_snap(system)
        ase_write(outfile, asesys, format=format)