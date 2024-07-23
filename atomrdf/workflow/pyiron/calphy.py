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
    method_dict["temperature"] = job.input.temperature
    method_dict["pressure"] = job.input.pressure

def add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["uri"] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software1 = {
        "uri": "http://demo.fiz-karlsruhe.de/matwerk/E447986",
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
            "label": "Pressure",
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
    outputs.append(
        {
            "label": "AtomicVolume",
            "value": np.round(job['output/atomic_volume'], decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }

    )    

    structure = job.get_structure(frame=-1)
    lx = np.linalg.norm(structure.cell[0])
    ly = np.linalg.norm(structure.cell[1])
    lz = np.linalg.norm(structure.cell[2])

    outputs.append(
        {
            "label": "SimulationCellLength_x",
            "value": np.round(lx, decimals=4),
            "unit": "ANGSTROM",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "SimulationCellLength_y",
            "value": np.round(ly, decimals=4),
            "unit": "ANGSTROM",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "SimulationCellLength_z",
            "value": np.round(lz, decimals=4),
            "unit": "ANGSTROM",
            "associate_to_sample": True,
        }
    )    
    method_dict['outputs'] =  outputs