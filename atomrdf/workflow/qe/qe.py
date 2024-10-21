"""
Wrappers for pyiron jobs
"""
import os
import numpy as np
import ast

from ase.io import read
from atomrdf.structure import System
from ase.io.espresso import read_fortran_namelist
from atomrdf.io import _convert_tab_to_dict

def _parse_inp(file):
    sample = None
    comments = []
    with open(file, 'r') as fin:
        for line in fin:
            line = line.strip()
            line = line.split('!')
            if len(line) > 1:
                comments.append(line[-1].strip())
    for comment in comments:
        rcomment = comment.split(':')
        if rcomment[0] == 'sample':
            sample = comment
            return sample
    return sample



def inform_graph(pr, kg):
    pass

def process_job(job):
    if not len(job)==2:
        raise ValueError('Job must be a tuple with two items: (quantum_espresso_input_file, quantum_espresso_output_file)')
    infile = job[0]
    outfile = job[1]

    method_dict = {}
    method_dict['intermediate'] = False
    get_structures(job, method_dict)
    identify_method(job, method_dict)
    add_software(method_dict)
    extract_calculated_quantities(job, method_dict)
    method_dict['path'] = os.path.abspath(os.path.dirname(infile))
    return method_dict

def get_structures(job, method_dict):
    infile = job[0]
    outfile = job[1]

    initial_ase_structure = read(infile, format='espresso-in')
    initial_pyscal_structure = System.read.ase(initial_ase_structure)

    #try to get initial sample id
    initial_sample_id = _parse_inp(infile)

    final_ase_structure = read(outfile, format='espresso-out')
    final_pyscal_structure = System.read.ase(final_ase_structure)
    method_dict['structure'] = {'initial': initial_pyscal_structure, 
                'final': final_pyscal_structure,} 
    method_dict['sample'] =  {'initial':initial_sample_id, 
                'final': None}

def identify_method(job, method_dict):
    infile = job[0]
    outfile = job[1]

    with open(infile, 'r') as fin:
        data, tab = read_fortran_namelist(fin)
    
    tab = _convert_tab_to_dict(tab)
    calc_method = data['control']['calculation']

    dof = []
    if calc_method in ['scf', 'nscf']:
        pass
    elif calc_method == 'relax':
        dof.append('AtomicPositionRelaxation')
    elif calc_method == 'vc-relax':
        dof.append('AtomicPositionRelaxation')
        dof.append('CellShapeRelaxation')
        dof.append('CellVolumeRelaxation')
    else:
        raise ValueError('Unknown calculation method')
    
    method = 'DensityFunctionalTheory'
    method_dict['method'] = method
    method_dict['dof'] = dof

    encut = data['system']['ecutwfc']
    #convert to eV
    inputs = []
    inputs.append(
        {
            "label": "EnergyCutoff",
            "value": encut*13.6057039763,
            "unit": "EV",
        }
    )

    #get kpoints
    if tab['K_POINTS']['extra'] == 'automatic':
        inputs.append(
            {
                "label": "MonkhorstPackKPointMesh",
                "value": " ".join(tab['K_POINTS']['value'][0].split()[:3]),
            }
        )
    elif tab['K_POINTS']['extra'] == 'gamma':
        inputs.append(
            {
                "label": "GammaCenteredKPointMesh",
                "value": " ",
            }
        )
    else:
        inputs.append(
            {
                "label": "ExplicitKPointMesh",
                "value": " ",
            }
        )
    method_dict['inputs'] = inputs

    #get pseudopotentials
    pseudo = None
    with open(outfile, 'r') as fin:
        for line in fin:
            if 'Exchange-correlation' in line:
                pseudo = line.split('=')[0].strip()
                break
    
    if pseudo is not None:
        method_dict['xc_functional'] = pseudo

def add_software(method_dict):
    software = {
        "uri": "https://www.quantum-espresso.org/",
        "label": "QuantumEspresso",
    }
    method_dict["software"] = [software]
    
def extract_calculated_quantities(job, method_dict):
    infile = job[0]
    outfile = job[1]   
    
    struct = read(outfile, format='espresso-out')
    outputs = []
    outputs.append(
        {
            "label": "TotalEnergy",
            "value": np.round(struct.get_total_energy(), decimals=5),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "Volume",
            "value": np.round(struct.get_volume(), decimals=5),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }
    ) 
