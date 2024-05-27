import os
import numpy as np
import ast
from atomrdf.structure import System

def process_job(job):
    method_dict = {}
    method_dict['intermediate'] = False
    get_structures(job, method_dict)
    identify_method(job, method_dict)
    extract_calculated_quantities(job, method_dict)
    add_software(method_dict)
    method_dict['path'] = get_simulation_folder(job)
    return method_dict

def get_simulation_folder(job):
    return os.path.join(job.project.path, f'{job.name}_hdf5')

def get_simulation_raw_folder(job):
    return os.path.join(job.project.path, f'{job.name}_hdf5', f'{job.name}')

def get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)
    initial_pyscal_structure = System.read.ase(initial_pyiron_structure)

    initial_sample_id = None
    
    if "sample_id" in initial_pyiron_structure.info.keys():
        initial_sample_id = initial_pyiron_structure.info["sample_id"]
    
    #now we can try to parse the POSCAR file directly here
    if initial_sample_id is None:
        #try to parse job directly, see if we have the structure written down in comments
        job.decompress()
        poscar_file_locations = [os.path.join(get_simulation_folder(job), 'POSCAR'),
                                os.path.join(get_simulation_raw_folder(job), 'POSCAR')]
        for poscar_file in poscar_file_locations:
            if os.path.exists(poscar_file):
                lines = []
                with open(poscar_file, 'r') as f:
                    for line in f:
                        lines.append(line)
                        break
                if 'sample' in lines[0]:
                    initial_sample_id = lines[0].strip()
                    break
                
    #add final structure
    final_pyscal_structure = System.read.ase(final_pyiron_structure)

    # now we do rthe transfer
    method_dict['structure'] = {'initial': initial_pyscal_structure, 
                'final': final_pyscal_structure,} 
    method_dict['sample'] =  {'initial':initial_sample_id, 
                'final': None}
    
def identify_method(job, method_dict):
    #get dof
    indf = job.input.incar.to_dict()
    params = indf['data_dict']['Parameter']
    vals = indf['data_dict']['Value']
    mlist = []
    for p,v in zip(params, vals):
        mlist.append(p + '=' + v)
    mstring = ';'.join(mlist)
    raw = mstring.split(';')
    mdict = {}
    for r in raw:
        rsplit = r.split('=')
        if len(rsplit) == 2:
            mdict[rsplit[0].replace(' ','')] = rsplit[1].replace(' ','')
    dof = []
    if 'ISIF' in mdict.keys():
        if mdict['ISIF'] in ['0', '1', '2']:
            dof.append('AtomicPositionRelaxation')
        elif mdict['ISIF'] == '3':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellShapeRelaxation')
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '4':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellShapeRelaxation')
        elif mdict['ISIF'] == '5':
            dof.append('CellShapeRelaxation')
        elif mdict['ISIF'] == '6':
            dof.append('CellShapeRelaxation')
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '7':
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '8':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellVolumeRelaxation')
    if 'NSW' in mdict.keys():
        if mdict['NSW'] == '0':
            dof = []

    method = 'DensityFunctionalTheory'
    method_dict['method'] = method
    method_dict['dof'] = dof

    encut = mdict['ENCUT'] 
    method_dict['encut'] = encut

    indf = job.input.to_dict()['kpoints/data_dict']
    params = indf['Parameter']
    vals = indf['Value']   

    kpoint_type = vals[2]
    kpoint_grid = vals[3]
    method_dict['kpoint_type'] = kpoint_type
    method_dict['kpoint_grid'] = kpoint_grid

    indf = job.input.to_dict()['potcar/data_dict']
    xc = indf['Value'][0]
    method_dict['xc_functional'] = xc

def add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["uri"] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software = {
        "uri": "https://www.vasp.at/",
        "label": "VASP",
    }
    method_dict["software"] = [software]

def extract_calculated_quantities(job, method_dict):
    """
    Extracts calculated quantities from a job.

    Parameters
    ----------
    job : pyiron.Job
        The job object containing the calculated quantities.

    Returns
    -------
    list
        A list of dictionaries, each containing the label, value, unit, and associate_to_sample of a calculated quantity.

    """
    outputs = []
    outputs.append(
        {
            "label": "TotalEnergy",
            "value": np.round(job.output.energy_tot[-1], decimals=5),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "TotalVolume",
            "value": np.round(job.output.volume[-1], decimals=5),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }
    )
    method_dict['outputs'] =  outputs