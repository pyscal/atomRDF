"""
Wrappers for pyiron jobs
"""

import os
from functools import partial, update_wrapper
from pyscal3.core import structure_dict, element_dict

from atomrdf.structure import _make_crystal
import atomrdf.workflow.pyiron.lammps as lammps

def process_job(job):
    """
    Checkes if the job is valid and creates the necessary output dict
    for the job.

    Parameters
    ----------
    job : pyiron.Job
        The pyiron job object to check.
    
    Raises
    ------
    TypeError
        If the job is not a valid pyiron job.
    """
    valid_jobs = [
        "Lammps",
    ]

    if type(job).__name__ == 'Lammps':
        return lammps.process_job(job)
    elif type(job).__name__ == 'Murnaghan':
        return process_murnaghan_job(job)
    else:
        raise TypeError("These type of pyiron Job is not currently supported")
    
    

def inform_graph(pr, kg):
    """
    this function in general can be used to do extra methods to set up things as needed
    for the workflow environment. 

    For example, for pyiron, this updates the project object to have the graph and creator objects
    """

    try:
        from pyiron_base import Creator, PyironFactory
        from pyiron_atomistics.atomistics.structure.atoms import (
            ase_to_pyiron,
            pyiron_to_ase,
        )
        import pyiron_atomistics.atomistics.structure.factory as sf
    except ImportError:
        raise ImportError("Please install pyiron_base and pyiron_atomistics")

    class AnnotatedStructureFactory:
        def __init__(self, graph):
            self._graph = graph

        def bulk(
            self,
            element,
            repetitions=None,
            crystalstructure=None,
            a=None,
            covera=None,
            cubic=True,
            graph=None,
            label=None,
        ):

            if crystalstructure is None:
                crystalstructure = element_dict[element]["structure"]
                if a is None:
                    a = element_dict[element]["lattice_constant"]

            struct = _make_crystal(
                crystalstructure,
                repetitions=repetitions,
                lattice_constant=a,
                ca_ratio=covera,
                element=element,
                primitive=not cubic,
                graph=self._graph,
                label=label,
            )

            ase_structure = struct.write.ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
            pyiron_structure.info["sample_id"] = struct.sample
            return pyiron_structure

        def grain_boundary(
            self,
            element,
            axis,
            sigma,
            gb_plane,
            repetitions=(1, 1, 1),
            crystalstructure=None,
            a=1,
            overlap=0.0,
            graph=None,
            label=None,
        ):

            struct = self._graph._annotated_make_grain_boundary(
                axis,
                sigma,
                gb_plane,
                structure=crystalstructure,
                element=element,
                lattice_constant=a,
                repetitions=repetitions,
                overlap=overlap,
                graph=self._graph,
                label=label,
            )

            ase_structure = struct.write.ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
            pyiron_structure.info["sample_id"] = struct.sample
            return pyiron_structure
        
    class StructureFactory(sf.StructureFactory):
        def __init__(self, graph):
            super().__init__()
            self._annotated_structure = AnnotatedStructureFactory(graph)

        @property
        def annotated_structure(self):
            return self._annotated_structure

    class StructureCreator(Creator):
        def __init__(self, project):
            super().__init__(project)
            self._structure = StructureFactory(project.graph)

        @property
        def structure(self):
            return self._structure

    pr.graph = kg
    pr._creator = StructureCreator(pr)

def process_murnaghan_job(job):
    #murnaghan job processing; add individual lammps jobs first
    job_dicts = []
    for jobid in job.child_ids:
        child_job = job.project.load(jobid)
        if type(child_job).__name__ == 'Lammps':
            single_job_dict = process_lammps_job(child_job)
            #note that we turn the jobs to child here
            single_job_dict['intermediate'] = True
            job_dicts.append(single_job_dict)
    #create an additional jobdict with the murnaghan job
    murnaghan_dict = {}
    murnaghan_dict['id'] = job.id
    murnaghan_dict['structure'] = get_structures(job)['structure']
    murnaghan_dict['sample'] = get_structures(job)['sample']
    murnaghan_dict['intermediate'] = False
    murnaghan_dict['path'] = get_simulation_folder(job)

    #add the murnaghan method
    murnaghan_dict['method'] = "EquationOfState"
    outputs = []
    outputs.append(
        {
            "label": "EquilibriumEnergy",
            "value": np.round(job['output/equilibrium_energy'], decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "EquilibriumVolume",
            "value": np.round(job['output/equilibrium_volume'], decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "BulkModulus",
            "value": np.round(job['output/equilibrium_bulk_modulus'], decimals=2),
            "unit": "GigaPA",
            "associate_to_sample": True,
        }
    )
    murnaghan_dict['outputs'] = outputs
    murnaghan_dict = add_software_lammps(murnaghan_dict)
    job_dicts.append(murnaghan_dict)
    return job_dicts



def add_software_vasp(mdict):
    mdict["workflow_manager"] = {}
    mdict["workflow_manager"]["uri"] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    mdict["workflow_manager"]["label"] = "pyiron"
    # and finally code details

    software = {
        "uri": "https://www.vasp.at/",
        "label": "VASP",
    }
    mdict["software"] = [software]
    return mdict



def vasp_extract_calculated_quantities(job):
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
    aen = job.output.energy_tot[-1]
    avol = job.output.volume[-1]
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
    return outputs



