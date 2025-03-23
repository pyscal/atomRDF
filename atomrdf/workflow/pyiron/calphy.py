import os
import numpy as np
import ast
from atomrdf.structure import System
import atomrdf.workflow.pyiron.lammps as lammps

def process_job(job):
    method_dict = {}
    method_dict['intermediate'] = False
    lammps.get_structures(job, method_dict)
    
    identify_method(job, method_dict)
    extract_calculated_quantities(job, method_dict)
    add_software(method_dict)
    get_simulation_folder(job, method_dict)
    return method_dict

def get_simulation_folder(job, method_dict):
    method_dict['path'] = os.path.join(job.project.path, f'{job.name}_hdf5')


def identify_method(job, method_dict):
    pressure = job.input.pressure
    if pressure is None:
        iso = True
        fix_lattice = True
    elif np.isscalar(pressure):
        iso = True
        fix_lattice = False
    elif np.shape(pressure) == (1,):
        iso = True
        fix_lattice = False
    elif np.shape(pressure) == (2,):
        iso = True
        fix_lattice = False
    elif np.shape(pressure) == (1, 3):
        iso = False
        fix_lattice = False
    elif np.shape(pressure) == (2, 3):
        iso = False
        fix_lattice = False
    
    dof = []
    dof.append("AtomicPositionRelaxation")
    ensemble = 'IsothermalIsobaricEnsemble'

    if not fix_lattice:
        dof.append("CellVolumeRelaxation")
        ensemble = "CanonicalEnsemble"

    if not iso:
        dof.append("CellShapeRelaxation")

    method_dict["dof"] = dof
    method_dict["ensemble"] = ensemble

    #now potential
    ps = job.potential.Config.values[0][0].strip().split('pair_style ')[-1]
    name = job.potential.Name.values[0]
    potstr = job.potential.Citations.values[0]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]

    method_dict["potential"] = {}
    method_dict["potential"]["type"] = ps
    method_dict["potential"]["label"] = name
    if url is not None:
        method_dict["potential"]["uri"] = url
    else:
        method_dict["potential"]["uri"] = name
    method_dict['method'] = 'ThermodynamicIntegration'
    method_dict['inputs'] = []
    method_dict['inputs'].append(
        {
            "label": "Pressure",
            "value": job.input.pressure,
            "unit": "BAR",
        }
    )
    method_dict['inputs'].append(
        {
            "label": "Temperature",
            "value": job.input.temperature,
            "unit": "K",
        }
    )    

def add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["uri"] = "https://doi.org/10.1016/j.commatsci.2018.07.043"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software1 = {
        "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
        "label": "LAMMPS",
    }

    software2 = {
        "uri": "https://doi.org/10.5281/zenodo.10527452",
        "label": "Calphy",
    }
    method_dict["software"] = [software1, software2]

def extract_calculated_quantities(job, method_dict):

    outputs = []
    outputs.append(
        {
            "label": "FreeEnergy",
            "value": np.round(job['output/energy_free'], decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "VirialPressure",
            "value": np.round(job['output/pressure'], decimals=4),
            "unit": "GigaPA",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "Temperature",
            "value": np.round(job['output/temperature'], decimals=2),
            "unit": "K",
            "associate_to_sample": True,
        }
    )  
    method_dict['outputs'] =  outputs