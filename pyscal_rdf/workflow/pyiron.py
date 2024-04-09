"""
Wrappers for pyiron jobs
"""
import os
from functools import partial, update_wrapper
import pyscal_rdf.workflow.workflow as wf
from pyscal_rdf.structure import _make_crystal
from pyscal3.core import structure_dict, element_dict


def _check_if_job_is_valid(job):
    valid_jobs = ['Lammps', ]
    
    if not type(job).__name__ in valid_jobs:
        raise TypeError('These type of pyiron Job is not currently supported')


def update_project(pr):
    """
    Update project to add extra creator functions
    """

    try:
        from pyiron_base import Creator, PyironFactory
        from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase
    except ImportError:
        raise ImportError('Please install pyiron_base and pyiron_atomistics')

    class StructureFactory(PyironFactory):
        def __init__(self):
            pass
    
        def bulk(self, 
            element,
            repetitions=None, 
            crystalstructure=None,
            a=None,
            covera=None,
            cubic=True,
            graph=None):

            if crystalstructure is None:
                crystalstructure = element_dict[element]['structure']
                if a is None:
                    a = element_dict[element]['lattice_constant']
        
            struct = _make_crystal(crystalstructure,
                repetitions=repetitions,
                lattice_constant=a,
                ca_ratio = covera,
                element = element,
                primitive = not cubic,
                graph=graph,
                )
            
            ase_structure = struct.write.ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
            if graph is not None:
                pyiron_structure.info['sample_id'] = struct.sample
            return pyiron_structure


        def grain_boundary(self,
            element,
            axis,
            sigma,
            gb_plane,
            repetitions = (1,1,1),
            crystalstructure=None,
            a=1,
            overlap=0.0,
            graph=None,
            ):

            struct = self._graph._annotated_make_grain_boundary(axis,
                sigma,
                gb_plane,
                structure = crystalstructure,
                element=element,
                lattice_constant=a,
                repetitions=repetitions,
                overlap=overlap,
                graph=graph)

            ase_structure = struct.write.ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
            if graph is not None:
                pyiron_structure.info['sample_id'] = struct.sample
            return pyiron_structure


    class StructureCreator(Creator):
        def __init__(self, project):
            super().__init__(project)
            self._structure = StructureFactory()
    
        @property
        def annotated_structure(self):
            return self._structure
    
    pr._creator = StructureCreator(pr)