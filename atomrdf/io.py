import mendeleev
import numpy as np
from ase.io.espresso import read_fortran_namelist
import os
import warnings

def write_espresso(s, inputfile, copy_from=None, pseudo_files=None):
    data = None
    if copy_from is not None:
        if os.path.exists(copy_from):
            try:
                with open(copy_from, 'r') as fin:
                    data, _ = read_fortran_namelist(fin)
            except:
                warnings.warn(f'Error reading {copy_from}, a clean file will be written')
                copy=True
    
    if data is None:
        data = {}
        data['system'] = {}
        data['control'] = {}
    
    lines = {}
    lines['CELL_PARAMETERS angstrom'] = []

    for vec in s.box:
        lines['CELL_PARAMETERS angstrom'].append(' '.join([str(x) for x in vec])) 

    cds = s.direct_coordinates
    species = s.atoms.species

    unique_species = np.unique(species)
    if pseudo_files is not None:
        if not len(pseudo_files) == len(unique_species):
            raise ValueError('Number of pseudo files must match number of unique species')
        pseudo_dirs = [os.path.dirname(os.path.abspath(pseudo_file)) for pseudo_file in pseudo_files]
        if not len(np.unique(pseudo_dirs)) == 1:
            raise ValueError('All pseudo files must be in the same directory')
        data['control']['pseudo_dir'] = f"'{pseudo_dirs[0]}'"
    else:
        pseudo_files = ['None' for x in range(len(unique_species))]

    lines['ATOMIC_SPECIES'] = []

    for count, us in enumerate(unique_species):
        chem = mendeleev.element(us)
        lines['ATOMIC_SPECIES'].append(f'{us} {chem.atomic_weight} {os.path.basename(pseudo_files[count])}')


    lines['ATOMIC_POSITIONS crystal'] = []

    for cd, sp in zip(cds, species):
        lines['ATOMIC_POSITIONS crystal'].append(f'{sp} {cd[0]} {cd[1]} {cd[2]}')

    data['system']['ibrav'] = 0
    data['system']['nat'] = len(species)
    data['system']['ntyp'] = len(unique_species)

    with open(inputfile, 'w') as fout:
        if s.sample is not None:
            fout.write(f'! {s.sample.toPython()}\n\n')
        for key, val in data.items():
            fout.write(f'&{key.upper()}\n')
            for k, v in val.items():
                fout.write(f'   {k} = {v},\n')
            fout.write('/\n')
            fout.write('\n')

        for key, val in lines.items():
            fout.write(key)
            fout.write('\n')
            for v in val:
                fout.write(v)
                fout.write('\n')
            fout.write('\n')   