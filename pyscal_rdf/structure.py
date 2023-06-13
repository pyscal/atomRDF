import numpy as np
from pyscal.core import System
from pyscal.crystal_structures import structure_creator, elements, structures
from pyscal_rdf.graph import RDFGraph
from pyscal.grain_boundary import GrainBoundary

class StructureGraph(RDFGraph):
    def __init__(self, graph_file=None):
        super().__init__(graph_file=graph_file)
        self._element_dict = elements
        self._structure_dict = structures
        
    def create_element(self, element, repetitions=(1,1,1), 
                       noise=0, add_to_graph=True, names=False):
        """
        Create elements
        """
        if element in self._element_dict.keys():
            structure = self._element_dict[element]['structure']
            sys = structure_creator(structure,
                        repetitions=repetitions,
                        noise=noise,
                        lattice_constant=self._element_dict[element]['lattice_constant'],
                        element = element)
            if add_to_graph:
                self.add_structure_to_graph(sys, names=names)
            return sys
    
    def create_structure(self, structure, 
                         lattice_constant = 1.00, 
                         repetitions = None, ca_ratio = 1.633, 
                         noise = 0, element=None,
                         add_to_graph=True, names=False):
        if structure in self._structure_dict.keys():
            sys = structure_creator(structure,
                        repetitions=repetitions,
                        noise=noise,
                        lattice_constant=lattice_constant,
                        element = element,
                        )
            if add_to_graph:
                self.add_structure_to_graph(sys, names = names)
            return sys
    
    def read_structure(self, filename, format="lammps-dump",
                      add_to_graph=True, names=False):
        sys = System(filename, format=format)
        if add_to_graph:
            self.add_structure_to_graph(sys, names=names)
        return sys
    
    def create_grain_boundary(self, axis, 
                              sigma, gb_plane,
                              structure=None,
                              element=None, 
                              lattice_constant=1,
                              repetitions=(1,1,1),
                              overlap=0.0,
                              add_to_graph=True,
                              names=False):
        gb = GrainBoundary()
        gb.create_grain_boundary(axis=axis, sigma=sigma, 
                                 gb_plane=gb_plane)

        #use standard creation routine
        if structure is not None:
            sys = gb.populate_grain_boundary(structure, 
                                             repetitions = repetitions,
                                             lattice_parameter = lattice_constant,
                                             overlap=overlap)
        elif element is not None:
            sys = gb.populate_grain_boundary(element, 
                                             repetitions=repetitions,
                                             overlap=overlap)
        else:
            raise ValueError("Either structure or element should be provided")
            
        #mapping of the system can be done
        self.add_structure_to_graph(sys, names=names)
        gb_dict = {"GBPlane": " ".join(np.array(gb_plane).astype(str)),
                  "RotationAxis": axis,
                  "MisorientationAngle": gb.theta,
                  "GBType": gb.find_gb_character(),
                  "sigma": gb.sigma,
                  }
        self.add_gb(gb_dict)
        return sys
            

