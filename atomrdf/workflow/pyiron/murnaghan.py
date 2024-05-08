import os
import numpy as np
import ast
from atomrdf.structure import System
import atomrdf.workflow.pyiron.lammps as lammps

def process_job(job):
    #murnaghan job processing; add individual lammps jobs first
    job_dicts = []
    for jobid in job.child_ids:
        child_job = job.project.load(jobid)
        if type(child_job).__name__ == 'Lammps':
            single_job_dict = lammps.process_job(child_job)
            #note that we turn the jobs to child here
            single_job_dict['intermediate'] = True
            job_dicts.append(single_job_dict)
    
    #create an additional jobdict with the murnaghan job
    murnaghan_dict = {}
    lammps.get_structures(job, murnaghan_dict)
    murnaghan_dict['intermediate'] = False
    lammps.get_simulation_folder(job, murnaghan_dict)

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
    lammps.add_software(murnaghan_dict)
    job_dicts.append(murnaghan_dict)
    return job_dicts