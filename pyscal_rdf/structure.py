from pyscal.core import System
from pyscal.crystal_structures import structure_creator, elements, structures
from pyscal_rdf.graph import RDFGraph

class StructureGraph(RDFGraph):
    def __init__(self, graph_file=None):
        super().__init__(graph_file=graph_file)
        self._element_dict = elements
        self._structure_dict = structures
        
    def create_element(self, element, repetitions=(1,1,1), 
                       noise=0, add_to_graph=True):
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
                self.add_structure_to_graph(sys)
            return sys
    
    def create_structure(self, structure, 
                         lattice_constant = 1.00, 
                         repetitions = None, ca_ratio = 1.633, 
                         noise = 0, element=None,
                         add_to_graph=True):
        if structure in self._structure_dict.keys():
            sys = structure_creator(structure,
                        repetitions=repetitions,
                        noise=noise,
                        lattice_constant=lattice_constant,
                        element = element)
            if add_to_graph:
                self.add_structure_to_graph(sys)
            return sys
    
    def read_structure(self, filename, format="lammps-dump",
                      add_to_graph=True):
        sys = System(filename, format=format)
        if add_to_graph:
            self.add_structure_to_graph(sys)
        return sys

