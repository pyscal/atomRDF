import networkx as nx
import matplotlib.pyplot as plt


class Network:
    """
    Network representation of CMSO
    """
    def __init__(self):
        self.g = nx.DiGraph()
    
    def add(self, sub, pred, obj, dtype=None):
        self.g.add_node(sub)
        self.g.add_node(pred)
        self.g.add_node(obj, dtype=dtype)            
        self.g.add_edge(sub, pred)
        self.g.add_edge(pred, obj)
    
    def draw(self):
        nx.draw(self.g, with_labels=True, font_weight='bold')
        
    def get_shortest_path(self, source, destination):
        path = nx.shortest_path(self.g, source=source, target=destination)
        return path
    
class OntologyNetwork(Network):
    def __init__(self):
        super().__init__()
        self.add("Sample", "hasMaterial", "Material")
        self.add("Material", "hasComposition", "ChemicalComposition")
        self.add("ChemicalComposition", "hasElementRatio", "ratio", dtype="string")

        self.add("Sample", "hasSimulationCell", "SimulationCell")
        self.add("SimulationCell", "hasVolume", "volume", dtype="float")
        self.add("Sample", "hasNumberOfAtoms", "numberofatoms", dtype="integer")

        self.add("SimulationCell", "hasLength", "SimulationCellLength")
        self.add("SimulationCellLength", "hasLength_x", "simulationcelllengthx", dtype="float")
        self.add("SimulationCellLength", "hasLength_y", "simulationcelllengthy", dtype="float")
        self.add("SimulationCellLength", "hasLength_z", "simulationcelllengthz", dtype="float")

        self.add("SimulationCell", "hasVector", "SimulationCellVectorA")
        self.add("SimulationCellVectorA", "hasComponent_x", "simulationcellvectorax", dtype="float")
        self.add("SimulationCellVectorA", "hasComponent_y", "simulationcellvectoray", dtype="float")
        self.add("SimulationCellVectorA", "hasComponent_z", "simulationcellvectoraz", dtype="float")
        self.add("SimulationCell", "hasVector", "SimulationCellVectorB")
        self.add("SimulationCellVectorB", "hasComponent_x", "simulationcellvectorbx", dtype="float")
        self.add("SimulationCellVectorB", "hasComponent_y", "simulationcellvectorby", dtype="float")
        self.add("SimulationCellVectorB", "hasComponent_z", "simulationcellvectorbz", dtype="float")
        self.add("SimulationCell", "hasVector", "SimulationCellVectorC")
        self.add("SimulationCellVectorC", "hasComponent_x", "simulationcellvectorcx", dtype="float")
        self.add("SimulationCellVectorC", "hasComponent_y", "simulationcellvectorcy", dtype="float")
        self.add("SimulationCellVectorC", "hasComponent_z", "simulationcellvectorcz", dtype="float")

        self.add("SimulationCell", "hasAngle", "SimulationCellAngle")
        self.add("SimulationCellAngle", "hasAngle_alpha", "simulationcellanglealpha", dtype="float")
        self.add("SimulationCellAngle", "hasAngle_beta", "simulationcellanglebeta", dtype="float")
        self.add("SimulationCellAngle", "hasAngle_gamma", "simulationcellanglegamma", dtype="float")

        self.add("Material", "hasStructure", "CrystalStructure")
        self.add("CrystalStructure", "hasAltName", "crystalstructurename", dtype="string")
        self.add("CrystalStructure", "hasSpaceGroup", "SpaceGroup")
        self.add("SpaceGroup", "hasSpaceGroupSymbol", "spacegroupsymbol", dtype="string")
        self.add("SpaceGroup", "hasSpaceGroupNumber", "spacegroupnumber", dtype="integer")

        self.add("CrystalStructure", "hasUnitCell", "UnitCell")
        self.add("UnitCell", "hasLattice", "BravaisLattice")
        self.add("BravaisLattice", "hasLatticeSystem", "latticesystem", dtype="string")
        self.add("UnitCell", "hasLatticeParameter", "LatticeParameter")
        self.add("LatticeParameter", "hasLength_x", "latticeparameterx", dtype="float")
        self.add("LatticeParameter", "hasLength_y", "latticeparametery", dtype="float")
        self.add("LatticeParameter", "hasLength_z", "latticeparameterz", dtype="float")
        self.add("UnitCell", "hasAngle", "LatticeAngle")
        self.add("LatticeAngle", "hasAngle_alpha", "anglealpha", dtype="float")
        self.add("LatticeAngle", "hasAngle_beta", "anglebeta", dtype="float")
        self.add("LatticeAngle", "hasAngle_gamma", "anglegamma", dtype="float")        
        