"""
StructureGraph is the central object in pyscal_rdf which combines all the functionality
of :py:class:`pyscal_rdf.graph.RDFGraph` along with easy structural creation routines.
"""
import numpy as np
import copy

from pyscal_rdf.rdfsystem import System
from pyscal3.crystal_structures import structure_creator, elements, structures
from pyscal_rdf.graph import RDFGraph
from pyscal3.grain_boundary import GrainBoundary

class StructureGraph(RDFGraph):
    def __init__(self, graph_file=None, 
        store="Memory", 
        store_file=None,
        identifier="default_graph"):
        
        super().__init__(graph_file=graph_file, store=store, store_file=store_file, identifier=identifier)
        self._element_dict = elements
        self._structure_dict = structures

    def create_element(self, element, repetitions=(1,1,1), 
                       noise=0, add_to_graph=True, names=True):
        """
        Create a crystal structure of the given element

        Parameters
        ----------
        element : string
            chemical symbol of the element

        repetitions : list of ints of len 3, optional
            of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
            default `[1, 1, 1]`.

        noise : float, optional
            If provided add normally distributed noise with standard deviation `noise` to the atomic positions.
        
        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs

        Returns
        -------
        System: :py:class:`pyscal.core.System`

        Notes
        -----
        Python module Mendelev is used to get the lattice constant with which the system will be constructed.
        If it is an hexagonal lattice, the ideal c/a ration will be used.

        """
        if element in self._element_dict.keys():
            structure = self._element_dict[element]['structure']
            sys = System(source=structure_creator(structure,
                        repetitions=repetitions,
                        noise=noise,
                        lattice_constant=self._element_dict[element]['lattice_constant'],
                        element = element))
            if add_to_graph:
                self.add_structure_to_graph(sys, names=names)
                #sys.sample = self.sample
                #sys._atom_ids = copy.copy(self._atom_ids)
                #sys.graph = self
            return sys
    
    def create_structure(self, structure, 
                         lattice_constant = 1.00, 
                         repetitions = None, ca_ratio = 1.633, 
                         noise = 0, element=None,
                         add_to_graph=True, names=True):
        """
        Create a crystal structure and return it as a System object.

        Parameters
        ----------
        structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
            type of the crystal structure

        lattice_constant : float, optional
            lattice constant of the crystal structure, default 1

        repetitions : list of ints of len 3, optional
            of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
            default `[1, 1, 1]`.

        ca_ratio : float, optional
            ratio of c/a for hcp structures, default 1.633

        noise : float, optional
            If provided add normally distributed noise with standard deviation `noise` to the atomic positions.

        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs
        
        element : string, optional
            The chemical element

        Returns
        -------
        System: pyscal System

        """
        if structure in self._structure_dict.keys():
            sys = System(source=structure_creator(structure,
                        repetitions=repetitions,
                        noise=noise,
                        lattice_constant=lattice_constant,
                        element = element,
                        ))
            if add_to_graph:
                self.add_structure_to_graph(sys, names = names)
                #sys.sample = self.sample
                #sys._atom_ids = copy.copy(self._atom_ids)
                #sys.graph = self
            return sys
    
    def read_structure(self, filename, format="lammps-dump",
                      add_to_graph=True, names=True):
        """
        Read an input file and return it as a System object.

        Parameters
        ----------
        filename: string
            name of the input file

        format: string
            format of the input file

        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs

        Returns
        -------
        System: pyscal System
        system will be populated with given atoms and simulation box

        """
        sys = System(filename, format=format)
        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)
            #sys.sample = self.sample
            #sys._atom_ids = copy.copy(self._atom_ids)
            #sys.graph = self
        return sys
    
    def create_grain_boundary(self, axis, 
                              sigma, gb_plane,
                              structure=None,
                              element=None, 
                              lattice_constant=1,
                              repetitions=(1,1,1),
                              overlap=0.0,
                              add_to_graph=True,
                              names=True):
        """
        Create a grain boundary structure and return it as a System object.

        Parameters
        ----------
        axis: list of ints of length 3
            The grain boundary axis

        sigma: int
            sigma value of the grain boundary
        
        gb_plane: list of ints of length 3
            The grain boundary plane

        structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
            type of the crystal structure

        element : string, optional
            The chemical element

        lattice_constant : float, optional
            lattice constant of the crystal structure, default 1

        repetitions : list of ints of len 3, optional
            of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
            default `[1, 1, 1]`.

        overlap: float, optional
            overlap between the two grains

        add_to_graph: bool, optinal
            If False, the created structure will not be added to the graph

        names: bool, optional
            If True, names will be used as IDs            

        Returns
        -------
        System: pyscal System

        """
        gb = GrainBoundary()
        gb.create_grain_boundary(axis=axis, sigma=sigma, 
                                 gb_plane=gb_plane)

        #use standard creation routine
        if structure is not None:
            sys = System(source=gb.populate_grain_boundary(structure, 
                                             repetitions = repetitions,
                                             lattice_parameter = lattice_constant,
                                             overlap=overlap))
        elif element is not None:
            sys = System(source=gb.populate_grain_boundary(element, 
                                             repetitions=repetitions,
                                             overlap=overlap))
        else:
            raise ValueError("Either structure or element should be provided")
            
        #mapping of the system can be done
        self.add_structure_to_graph(sys, names=names)
        #sys.sample = self.sample
        #sys._atom_ids = copy.copy(self._atom_ids)
        #sys.graph = self
        gb_dict = {"GBPlane": " ".join(np.array(gb_plane).astype(str)),
                  "RotationAxis": axis,
                  "MisorientationAngle": gb.theta,
                  "GBType": gb.find_gb_character(),
                  "sigma": gb.sigma,
                  }
        self.add_gb(gb_dict)
        return sys
            

