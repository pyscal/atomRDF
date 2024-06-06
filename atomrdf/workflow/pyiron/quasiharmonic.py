import os
import numpy as np
import ast
from atomrdf.structure import System
import atomrdf.workflow.pyiron.lammps as lammps

def process_job(job):
    #murnaghan job processing; add individual lammps jobs first
    job_dicts = []
    
    #create an additional jobdict with the murnaghan job
    quasi_dict = {}
    lammps.get_structures(job, quasi_dict)
    quasi_dict['intermediate'] = False
    lammps.get_simulation_folder(job, quasi_dict)

    #add the murnaghan method
    quasi_dict['method'] = "QuasiHarmonicModel"
    outputs = []
    outputs.append(
        {
            "label": "QuasiHarmonicFreeEnergy",
            "value": np.round(job['output/free_energy'].T, decimals=4),
            "unit": "EV",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "QuasiHarmonicVolume",
            "value": np.round(job['output/volumes'].T, decimals=4),
            "unit": "ANGSTROM3",
            "associate_to_sample": True,
        }
    )
    outputs.append(
        {
            "label": "QuasiHarmonicTemperature",
            "value": np.round(job['output/temperatures'][0], decimals=2),
            "unit": "K",
            "associate_to_sample": True,
        }
    )
    quasi_dict['outputs'] = outputs
    lammps.add_software(quasi_dict)
    job_dicts.append(quasi_dict)
    return job_dicts