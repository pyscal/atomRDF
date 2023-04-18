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
        
    def get_shortest_path(self, source, target):
        if not source
        path = nx.shortest_path(self.g, source=source, target=target)
        return path
    
class OntologyNetwork(Network):
    def __init__(self):
        super().__init__()
        self.add("Sample", "hasMaterial", "Material")
        self.add("Material", "hasComposition", "ChemicalComposition")
        self.add("ChemicalComposition", "hasElementRatio", "ElementRatio", dtype="string")

        self.add("Sample", "hasSimulationCell", "SimulationCell")
        self.add("SimulationCell", "hasVolume", "Volume", dtype="float")
        self.add("Sample", "hasNumberOfAtoms", "NumberOfAtoms", dtype="integer")

        self.add("SimulationCell", "hasLength", "SimulationCellLength")
        self.add("SimulationCellLength", "hasLength_x", "SimulationCellLength_x", dtype="float")
        self.add("SimulationCellLength", "hasLength_y", "SimulationCellLength_y", dtype="float")
        self.add("SimulationCellLength", "hasLength_z", "SimulationCellLength_z", dtype="float")

        self.add("SimulationCell", "hasVector", "SimulationCellVectorA")
        self.add("SimulationCellVectorA", "hasComponent_x", "SimulationCellVectorA_x", dtype="float")
        self.add("SimulationCellVectorA", "hasComponent_y", "SimulationCellVectorA_y", dtype="float")
        self.add("SimulationCellVectorA", "hasComponent_z", "SimulationCellVectorA_z", dtype="float")
        self.add("SimulationCell", "hasVector", "SimulationCellVectorB")
        self.add("SimulationCellVectorB", "hasComponent_x", "SimulationCellVectorB_x", dtype="float")
        self.add("SimulationCellVectorB", "hasComponent_y", "SimulationCellVectorB_y", dtype="float")
        self.add("SimulationCellVectorB", "hasComponent_z", "SimulationCellVectorB_z", dtype="float")
        self.add("SimulationCell", "hasVector", "SimulationCellVectorC")
        self.add("SimulationCellVectorC", "hasComponent_x", "SimulationCellVectorC_x", dtype="float")
        self.add("SimulationCellVectorC", "hasComponent_y", "SimulationCellVectorC_y", dtype="float")
        self.add("SimulationCellVectorC", "hasComponent_z", "SimulationCellVectorC_z", dtype="float")

        self.add("SimulationCell", "hasAngle", "SimulationCellAngle")
        self.add("SimulationCellAngle", "hasAngle_alpha", "SimulationCellAngle_alpha", dtype="float")
        self.add("SimulationCellAngle", "hasAngle_beta", "SimulationCellAngle_beta", dtype="float")
        self.add("SimulationCellAngle", "hasAngle_gamma", "SimulationCellAngle_gamma", dtype="float")

        self.add("Material", "hasStructure", "CrystalStructure")
        self.add("CrystalStructure", "hasAltName", "CrystalStructureAltName", dtype="string")
        self.add("CrystalStructure", "hasSpaceGroup", "SpaceGroup")
        self.add("SpaceGroup", "hasSpaceGroupSymbol", "SpaceGroupSymbol", dtype="string")
        self.add("SpaceGroup", "hasSpaceGroupNumber", "SpaceGroupNumber", dtype="integer")

        self.add("CrystalStructure", "hasUnitCell", "UnitCell")
        self.add("UnitCell", "hasLattice", "BravaisLattice")
        self.add("BravaisLattice", "hasLatticeSystem", "LatticeSystem", dtype="string")
        self.add("UnitCell", "hasLatticeParameter", "LatticeParameter")
        self.add("LatticeParameter", "hasLength_x", "LatticeParameter_x", dtype="float")
        self.add("LatticeParameter", "hasLength_y", "LatticeParameter_y", dtype="float")
        self.add("LatticeParameter", "hasLength_z", "LatticeParameter_z", dtype="float")
        self.add("UnitCell", "hasAngle", "LatticeAngle")
        self.add("LatticeAngle", "hasAngle_alpha", "LatticeAngle_alpha", dtype="float")
        self.add("LatticeAngle", "hasAngle_beta", "LatticeAngle_beta", dtype="float")
        self.add("LatticeAngle", "hasAngle_gamma", "LatticeAngle_gamma", dtype="float")
        
    def get_path_from_sample(self, target):
        path = self.get_shortest_path(source="Sample", target=target)
        triplets = []
        for x in range(len(path)//2):
            triplets.append(path[2*x:2*x+3])
        return triplets
        
        