"""
Wrappers for pyiron jobs
"""
import os
import numpy as np
import ast

from ase.io import read
from atomrdf.structure import System
from ase.io.espresso import read_fortran_namelist

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
        data, _ = read_fortran_namelist(fin)
    
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
    method_dict['encut'] = encut*13.6057039763

    
