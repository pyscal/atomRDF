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
        self.g.add("Sample", "hasMaterial", "Material")
        self.g.add("Material", "hasComposition", "ChemicalComposition")
        self.g.add("ChemicalComposition", "hasElementRatio", "ratio", dtype="string")

        self.g.add("Sample", "hasSimulationCell", "SimulationCell")
        self.g.add("SimulationCell", "hasVolume", "volume", dtype="float")
        self.g.add("Sample", "hasNumberOfAtoms", "numberofatoms", dtype="integer")

        self.g.add("SimulationCell", "hasLength", "SimulationCellLength")
        self.g.add("SimulationCellLength", "hasLength_x", "simulationcelllengthx", dtype="float")
        self.g.add("SimulationCellLength", "hasLength_y", "simulationcelllengthy", dtype="float")
        self.g.add("SimulationCellLength", "hasLength_z", "simulationcelllengthz", dtype="float")

        self.g.add("SimulationCell", "hasVector", "SimulationCellVectorA")
        self.g.add("SimulationCellVectorA", "hasComponent_x", "simulationcellvectorax", dtype="float")
        self.g.add("SimulationCellVectorA", "hasComponent_y", "simulationcellvectoray", dtype="float")
        self.g.add("SimulationCellVectorA", "hasComponent_z", "simulationcellvectoraz", dtype="float")
        self.g.add("SimulationCell", "hasVector", "SimulationCellVectorB")
        self.g.add("SimulationCellVectorB", "hasComponent_x", "simulationcellvectorbx", dtype="float")
        self.g.add("SimulationCellVectorB", "hasComponent_y", "simulationcellvectorby", dtype="float")
        self.g.add("SimulationCellVectorB", "hasComponent_z", "simulationcellvectorbz", dtype="float")
        self.g.add("SimulationCell", "hasVector", "SimulationCellVectorC")
        self.g.add("SimulationCellVectorC", "hasComponent_x", "simulationcellvectorcx", dtype="float")
        self.g.add("SimulationCellVectorC", "hasComponent_y", "simulationcellvectorcy", dtype="float")
        self.g.add("SimulationCellVectorC", "hasComponent_z", "simulationcellvectorcz", dtype="float")

        self.g.add("SimulationCell", "hasAngle", "SimulationCellAngle")
        self.g.add("SimulationCellAngle", "hasAngle_alpha", "simulationcellanglealpha", dtype="float")
        self.g.add("SimulationCellAngle", "hasAngle_beta", "simulationcellanglebeta", dtype="float")
        self.g.add("SimulationCellAngle", "hasAngle_gamma", "simulationcellanglegamma", dtype="float")

        self.g.add("Material", "hasStructure", "CrystalStructure")
        self.g.add("CrystalStructure", "hasAltName", "crystalstructurename", dtype="string")
        self.g.add("CrystalStructure", "hasSpaceGroup", "SpaceGroup")
        self.g.add("SpaceGroup", "hasSpaceGroupSymbol", "spacegroupsymbol", dtype="string")
        self.g.add("SpaceGroup", "hasSpaceGroupNumber", "spacegroupnumber", dtype="integer")

        self.g.add("CrystalStructure", "hasUnitCell", "UnitCell")
        self.g.add("UnitCell", "hasLattice", "BravaisLattice")
        self.g.add("BravaisLattice", "hasLatticeSystem", "latticesystem", dtype="string")
        self.g.add("UnitCell", "hasLatticeParameter", "LatticeParameter")
        self.g.add("LatticeParameter", "hasLength_x", "latticeparameterx", dtype="float")
        self.g.add("LatticeParameter", "hasLength_y", "latticeparametery", dtype="float")
        self.g.add("LatticeParameter", "hasLength_z", "latticeparameterz", dtype="float")
        self.g.add("UnitCell", "hasAngle", "LatticeAngle")
        self.g.add("LatticeAngle", "hasAngle_alpha", "anglealpha", dtype="float")
        self.g.add("LatticeAngle", "hasAngle_beta", "anglebeta", dtype="float")
        self.g.add("LatticeAngle", "hasAngle_gamma", "anglegamma", dtype="float")        
        