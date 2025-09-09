import os
import numpy as np
import ast
from atomrdf.datamodels.workflow.workflow import Simulation
from atomrdf.datamodels.workflow.property import InputParameter, OutputParameter, CalculatedProperty

def ase_to_pyiron(structure):
    try:
        from pyiron_atomistics.atomistics.structure.atoms import (
            ase_to_pyiron,
            pyiron_to_ase,
        )
    except ImportError:
        raise ImportError("Please install pyiron_atomistics")
    
    pyiron_structure = ase_to_pyiron(structure)
    if 'id' in structure.info:
        pyiron_structure.info['id'] = structure.info['id']
    return pyiron_structure

def extract(job):
    if type(job).__name__ == 'Lammps':
        return process_job_lammps(job)

def process_job_lammps(job):
    method_dict = Simulation.template()
    identify_method(job, method_dict)
    add_software(method_dict)
    get_simulation_folder(job, method_dict)
    extract_calculated_quantities(job, method_dict)
    return method_dict

def get_simulation_folder(job, method_dict):
    method_dict['path'] = os.path.join(job.project.path, f'{job.name}_hdf5')

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

    method_dict["method"] = {'basename': md_method}
    method_dict["degrees_of_freedom"] = dof
    method_dict["thermodynamic_ensemble"] = {'basename': ensemble}

    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]

    pssplit = ps.split('/')
    if len(pssplit) > 1:
        ps = pssplit[0]

    method_dict["interatomic_potential"] = {
        'potential_type': ps,
    }

    if url is not None:
        method_dict["interatomic_potential"]["uri"] = url
    else:
        method_dict["interatomic_potential"]["uri"] = name

    #add temperature and pressure as inputs
    if temp is not None:
        method_dict['input_parameter'].append(
            {
                "basename": "Temperature",
                "label": "Temperature",
                "value": temp,
                "unit": "K",
            }
        )
    if press is not None:
        method_dict['input_parameter'].append(
            {
                "basename": "Pressure",
                "label": "Pressure",
                "value": press,
                "unit": "GigaPA",
            }
        )

def add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["uri"] = "https://doi.org/10.1016/j.commatsci.2018.07.043"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software = {
        "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
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
    energy_tot = np.mean(job.output.energy_tot)
    energy_pot = np.mean(job.output.energy_pot)
    energy_kin = energy_tot - energy_pot

    volume = np.mean(job.output.volume)
    
    outputs = []
    outputs.append(
        {
            "label": "TotalEnergy",
            "basename": "TotalEnergy",
            "value": np.round(energy_tot, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "PotentialEnergy",
            "basename": "PotentialEnergy",
            "value": np.round(energy_pot, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "KineticEnergy",
            "basename": "KineticEnergy",
            "value": np.round(energy_kin, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "SimulationCellVolume",
            "value": np.round(volume, decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
            "basename": "SimulationCellVolume",
        }
    )
    
    method_dict['calculated_property'] =  outputs