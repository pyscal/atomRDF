import os
import numpy as np
import ast
from atomrdf.datamodels.workflow.workflow import Simulation
from atomrdf.datamodels.workflow.property import (
    InputParameter,
    OutputParameter,
    CalculatedProperty,
)


def ase_to_pyiron(structure):
    try:
        from pyiron_atomistics.atomistics.structure.atoms import (
            ase_to_pyiron,
            pyiron_to_ase,
        )
    except ImportError:
        raise ImportError("Please install pyiron_atomistics")

    pyiron_structure = ase_to_pyiron(structure)
    if "id" in structure.info:
        pyiron_structure.info["id"] = structure.info["id"]
    return pyiron_structure


def extract(job):
    if type(job).__name__ == "Lammps":
        return lammps_process_job(job)
    elif type(job).__name__ == "Calphy":
        return calphy_process_job(job)
    elif type(job).__name__ == "Murnaghan":
        return murnaghan_process_job(job)
    elif type(job).__name__ == "QuasiHarmonicJob":
        return quasiharmonic_process_job(job)
    else:
        raise TypeError("These type of pyiron Job is not currently supported")


def lammps_process_job(job):
    method_dict = Simulation.template()
    _lammps_identify_method(job, method_dict)
    _lammps_add_software(method_dict)
    _lammps_get_simulation_folder(job, method_dict)
    _lammps_extract_calculated_quantities(job, method_dict)
    _lammps_get_structures(job, method_dict)
    return method_dict


def calphy_process_job(job):
    method_dict = Simulation.template()
    _calphy_identify_method(job, method_dict)
    _calphy_add_software(method_dict)
    _lammps_get_simulation_folder(job, method_dict)
    _calphy_extract_calculated_quantities(job, method_dict)
    _lammps_get_structures(job, method_dict)
    return method_dict


def murnaghan_process_job(job):
    # murnaghan job processing; add individual lammps jobs first
    job_dicts = []
    for jobid in job.child_ids:
        child_job = job.project.load(jobid)
        if type(child_job).__name__ == "Lammps":
            single_job_dict = lammps_process_job(child_job)
            # note that we turn the jobs to child here
            job_dicts.append(single_job_dict)

    # create an additional jobdict with the murnaghan job
    murnaghan_dict = Simulation.template()
    _lammps_get_structures(job, murnaghan_dict)
    _lammps_get_simulation_folder(job, murnaghan_dict)

    # add the murnaghan method
    murnaghan_dict["algorithm"] = "EquationOfStateFit"
    murnaghan_dict["method"] = "MolecularStatics"
    outputs = []
    outputs.append(
        {
            "label": "EquilibriumEnergy",
            "value": np.round(job["output/equilibrium_energy"], decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
            "basename": "TotalEnergy",
        }
    )
    outputs.append(
        {
            "label": "EquilibriumVolume",
            "value": np.round(job["output/equilibrium_volume"], decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
            "basename": "Volume",
        }
    )
    outputs.append(
        {
            "label": "TotalEnergy",
            "value": np.round(job["output/energy"], decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
            "basename": "TotalEnergy",
        }
    )
    outputs.append(
        {
            "label": "SimulationCellVolume",
            "value": np.round(job["output/volume"], decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
            "basename": "SimulationCellVolume",
        }
    )
    outputs.append(
        {
            "label": "BulkModulus",
            "basename": "BulkModulus",
            "value": np.round(job["output/equilibrium_bulk_modulus"], decimals=2),
            "unit": "GigaPA",
            "associate_to_sample": True,
        }
    )

    murnaghan_dict["calculated_property"] = outputs
    _lammps_add_software(murnaghan_dict)
    job_dicts.append(murnaghan_dict)
    return job_dicts


def quasiharmonic_process_job(job):
    # murnaghan job processing; add individual lammps jobs first
    job_dicts = []

    # create an additional jobdict with the murnaghan job
    quasi_dict = Simulation.template()
    _lammps_get_structures(job, quasi_dict)
    _lammps_get_simulation_folder(job, quasi_dict)

    # add the murnaghan method
    quasi_dict["algorithm"] = "QuasiHarmonicApproximation"
    quasi_dict["method"] = "MolecularStatics"
    outputs = []
    outputs.append(
        {
            "label": "QuasiharmonicFreeEnergy",
            "value": np.round(job["output/free_energy"].T, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
            "basename": "FreeEnergy",
        }
    )
    outputs.append(
        {
            "label": "QuasiharmonicVolume",
            "value": np.round(job["output/volumes"].T, decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
            "basename": "SimulationCellVolume",
        }
    )
    outputs.append(
        {
            "label": "QuasiharmonicTemperature",
            "value": np.round(job["output/temperatures"][0], decimals=2),
            "unit": "K",
            "associate_to_sample": True,
            "basename": "Temperature",
        }
    )
    quasi_dict["calculated_property"] = outputs
    _lammps_add_software(quasi_dict)
    job_dicts.append(quasi_dict)
    return job_dicts


def _lammps_get_simulation_folder(job, method_dict):
    method_dict["path"] = os.path.join(job.project.path, f"{job.name}_hdf5")


def _lammps_identify_method(job, method_dict):
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
        raw = input_dict["fix___ensemble"].split().format(i)
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

    method_dict["method"] = {"basename": md_method}
    method_dict["degrees_of_freedom"] = dof
    if ensemble is not None:
        method_dict["thermodynamic_ensemble"] = {"basename": ensemble}

    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]

    pssplit = ps.split("/")
    if len(pssplit) > 1:
        ps = pssplit[0]

    method_dict["interatomic_potential"] = {
        "potential_type": ps,
    }

    if url is not None:
        method_dict["interatomic_potential"]["uri"] = url
    else:
        method_dict["interatomic_potential"]["uri"] = name

    # add temperature and pressure as inputs
    if temp is not None:
        method_dict["input_parameter"].append(
            {
                "basename": "Temperature",
                "label": "Temperature",
                "value": temp,
                "unit": "K",
            }
        )
    if press is not None:
        method_dict["input_parameter"].append(
            {
                "basename": "Pressure",
                "label": "Pressure",
                "value": press,
                "unit": "GigaPA",
            }
        )


def _lammps_add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"][
        "uri"
    ] = "https://doi.org/10.1016/j.commatsci.2018.07.043"
    method_dict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software = {
        "uri": "https://doi.org/10.1016/j.cpc.2021.108171",
        "label": "LAMMPS",
    }
    method_dict["software"] = [software]


def _lammps_extract_calculated_quantities(job, method_dict):
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

    method_dict["calculated_property"] = outputs


def _lammps_get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)

    method_dict["initial_sample"] = initial_pyiron_structure
    method_dict["final_sample"] = final_pyiron_structure


def _calphy_identify_method(job, method_dict):
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
    ensemble = "IsothermalIsobaricEnsemble"

    if not fix_lattice:
        dof.append("CellVolumeRelaxation")
        ensemble = "CanonicalEnsemble"

    if not iso:
        dof.append("CellShapeRelaxation")

    method_dict["method"] = {"basename": "ThermodynamicIntegration"}
    method_dict["degrees_of_freedom"] = dof
    method_dict["thermodynamic_ensemble"] = {"basename": ensemble}

    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]

    pssplit = ps.split("/")
    if len(pssplit) > 1:
        ps = pssplit[0]

    method_dict["interatomic_potential"] = {
        "potential_type": ps,
    }

    if url is not None:
        method_dict["interatomic_potential"]["uri"] = url
    else:
        method_dict["interatomic_potential"]["uri"] = name

    # add temperature and pressure as inputs
    method_dict["input_parameter"].append(
        {
            "basename": "Temperature",
            "label": "Temperature",
            "value": job.input.temperature,
            "unit": "K",
        }
    )

    method_dict["input_parameter"].append(
        {
            "basename": "Pressure",
            "label": "Pressure",
            "value": job.input.pressure,
            "unit": "GigaPA",
        }
    )


def _calphy_add_software(method_dict):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"][
        "uri"
    ] = "https://doi.org/10.1016/j.commatsci.2018.07.043"
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


def _calphy_extract_calculated_quantities(job, method_dict):

    outputs = []
    outputs.append(
        {
            "label": "FreeEnergy",
            "basename": "FreeEnergy",
            "value": np.round(job["output/energy_free"], decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "VirialPressure",
            "basename": "VirialPressure",
            "value": np.round(job["output/pressure"], decimals=4),
            "unit": "GigaPA",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "Temperature",
            "basename": "Temperature",
            "value": np.round(job["output/temperature"], decimals=2),
            "unit": "K",
            "associate_to_sample": True,
        }
    )
    method_dict["calculated_property"] = outputs
