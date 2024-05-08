"""
LAMMPS specific functions for parsing

Use this a reference for specific implementations
"""
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
    get_simulation_folder(job, method_dict)
    return method_dict

def get_simulation_folder(job, method_dict):
    method_dict['path'] = os.path.join(job.project.path, f'{job.name}_hdf5')

def get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)
    initial_pyscal_structure = System.read.ase(initial_pyiron_structure)
    initial_sample_id = None

    if "sample_id" in initial_pyiron_structure.info.keys():
        initial_sample_id = initial_pyiron_structure.info["sample_id"]

    #add final structure
    final_pyscal_structure = System.read.ase(final_pyiron_structure)

    # now we do rthe transfer
    method_dict['structure'] = {'initial': initial_pyscal_structure, 
                'final': final_pyscal_structure,} 
    method_dict['sample'] =  {'initial':initial_sample_id, 
                'final': None}

def identify_method(job, method_dict):
    job_dict = job.input.to_dict()
    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
    }
    dof = []
    temp = None
    press = None
    md_method = None
    ensemble = None

    if "min_style" in input_dict.keys():
        dof.append("AtomicPositionRelaxation")
        dof.append("CellVolumeRelaxation")
        md_method = "MolecularStatics"

    elif "nve" in input_dict["fix___ensemble"]:
        if int(input_dict["run"]) == 0:
            method = "static"
            md_method = "MolecularStatics"
            ensemble = "MicrocanonicalEnsemble"

        elif int(input_dict["run"]) > 0:
            method = "md_nve"
            dof.append("AtomicPositionRelaxation")
            md_method = "MolecularDynamics"
            ensemble = "MicrocanonicalEnsemble"

    elif "nvt" in input_dict["fix___ensemble"]:
        method = "md_nvt"
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        dof.append("AtomicPositionRelaxation")
        md_method = "MolecularDynamics"
        ensemble = "CanonicalEnsemble"

    elif "npt" in input_dict["fix___ensemble"]:
        dof.append("AtomicPositionRelaxation")
        dof.append("CellVolumeRelaxation")
        if "aniso" in input_dict["fix___ensemble"]:
            method = "md_npt_aniso"
            dof.append("CellShapeRelaxation")
        else:
            method = "md_npt_iso"
        md_method = "MolecularDynamics"
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        press = float(raw[7])
        ensemble = "IsothermalIsobaricEnsemble"

    method_dict["method"] = md_method
    method_dict["temperature"] = temp
    method_dict["pressure"] = press
    method_dict["dof"] = dof
    method_dict["ensemble"] = ensemble

    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
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

def add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["uri"] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software = {
        "uri": "http://demo.fiz-karlsruhe.de/matwerk/E447986",
        "label": "LAMMPS",
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
    aen = np.mean(job.output.energy_tot)
    avol = np.mean(job.output.volume)
    outputs = []
    outputs.append(
        {
            "label": "TotalEnergy",
            "value": np.round(aen, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "TotalVolume",
            "value": np.round(avol, decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }
    )
    method_dict['outputs'] =  outputs