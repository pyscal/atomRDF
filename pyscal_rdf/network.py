import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class Network:
    """
    Network representation of CMSO
    """
    def __init__(self):
        self.g = nx.DiGraph()
    
    def add(self, sub, pred, obj, dtype=None, pred_prefix="cmso"):
        pred = f'{pred_prefix}:{pred}'
        self.g.add_node(sub, node_type="object")
        self.g.add_node(pred, node_type="property")
        if dtype is not None:
            nd = "data"
        else:
            nd = "object"
        self.g.add_node(obj, dtype=dtype, node_type=nd)            
        self.g.add_edge(sub, pred)
        self.g.add_edge(pred, obj)
    
    def draw(self):
        nx.draw(self.g, with_labels=True, font_weight='bold')
        
    def get_shortest_path(self, source, target):
        path = nx.shortest_path(self.g, source=source, target=target)
        return path
    
class OntologyNetwork(Network):
    def __init__(self):
        super().__init__()
        self.add("Sample", "hasMaterial", "Material")
        self.add("Material", "hasElementRatio", "ElementRatio", dtype="string")

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
        self.add("CrystalStructure", "hasSpaceGroupSymbol", "SpaceGroupSymbol", dtype="string")
        self.add("CrystalStructure", "hasSpaceGroupNumber", "SpaceGroupNumber", dtype="integer")

        self.add("CrystalStructure", "hasUnitCell", "UnitCell")
        self.add("UnitCell", "hasBravaisLattice", "LatticeSystem")
        self.add("UnitCell", "hasLatticeParameter", "LatticeParameter")
        self.add("LatticeParameter", "hasLength_x", "LatticeParameter_x", dtype="float")
        self.add("LatticeParameter", "hasLength_y", "LatticeParameter_y", dtype="float")
        self.add("LatticeParameter", "hasLength_z", "LatticeParameter_z", dtype="float")
        self.add("UnitCell", "hasAngle", "LatticeAngle")
        self.add("LatticeAngle", "hasAngle_alpha", "LatticeAngle_alpha", dtype="float")
        self.add("LatticeAngle", "hasAngle_beta", "LatticeAngle_beta", dtype="float")
        self.add("LatticeAngle", "hasAngle_gamma", "LatticeAngle_gamma", dtype="float")

        #add GB properties
        self.add("Material", "hasDefect", "Defect", pred_prefix="cmso")
        self.add("Defect", "type", "GrainBoundary", pred_prefix="rdf")
        self.add("Defect", "type", "TwistBoundary", pred_prefix="rdf")
        self.add("Defect", "type", "TiltBoundary", pred_prefix="rdf")
        self.add("Defect", "type", "SymmetricTiltBoundary", pred_prefix="rdf")
        self.add("Defect", "type", "MixedBoundary", pred_prefix="rdf")
        self.add("Defect", "hasSigmaValue", "Sigma", dtype="integer", pred_prefix="pldo")
        self.add("Defect", "hasGBPlane", "GBPlane", pred_prefix="pldo", dtype="string")
        self.add("Defect", "hasRotationAxis", "RotationAxis", pred_prefix="pldo", dtype="string")
        self.add("Defect", "hasMisorientationAngle", "MisorientationAngle", pred_prefix="pldo", dtype="float")

        #add vacancy
        self.add("Defect", "type", "Vacancy", pred_prefix="rdf")
        self.add("SimulationCell", "hasVacancyConcentration", "VacancyConcentration", pred_prefix="podo", dtype="float")
        self.add("SimulationCell", "hasNumberOfVacancies", "NumberOfVacancy", pred_prefix="podo", dtype="integer")
         

    def get_path_from_sample(self, target):
        path = self.get_shortest_path(source="Sample", target=target)
        triplets = []
        for x in range(len(path)//2):
            triplets.append(path[2*x:2*x+3])
        return triplets
        
    def formulate_query(self, target, value):
        #first get triplets
        triplets = self.get_path_from_sample(target)
        #start building query
        query = self._formulate_query_path(triplets)
        query.append(self._formulate_filter_expression(triplets, value))
        query.append("}")
        query = " ".join(query)
        return query
        
    
    def _formulate_query_path(self, triplets):
        query = []
        query.append("PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>")
        query.append("PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>")
        query.append("PREFIX podo: <https://purls.helmholtz-metadaten.de/podo/>")
        query.append("PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>")
        query.append("SELECT DISTINCT ?sample")
        query.append("WHERE {")
        for triple in triplets:
            query.append("    ?%s %s ?%s ."%(triple[0].lower(), 
                                                  triple[1], 
                                                  triple[2].lower()))
        return query
    
    def _formulate_filter_expression(self, triplets, value):                       
        value, datatype = self._check_value(value)      
        last_val = self.g.nodes[triplets[-1][-1]]
        last_val_name = triplets[-1][-1].lower()
        
        #if it is nodetype data
        if last_val['node_type'] == "data":
            if datatype == "multi_string":
                qstr = self._formulate_or_string_query(last_val, 
                                                   last_val_name, 
                                                   value)
            elif datatype == "multi_number":
                qstr = self._formulate_range_number_query(last_val, 
                                                   last_val_name, 
                                                   value)
            else:
                qstr = self._formulate_equal_query(last_val, 
                                                   last_val_name, 
                                                   value)
            return qstr
        else:
            raise NotImplementedError("Non-data queries are not implemented")
    
    def _check_value(self, value):
        if isinstance(value, list):
            if not len(value) == 2:
                raise ValueError("value can be maximum length 2")
        else:
            value = [value]
        if all(isinstance(x, str) for x in value):
            datatype = "string"
        elif all(isinstance(x, (int, float)) for x in value):
            datatype = "number"
        else:
            raise TypeError("Values have to be of same type")
        if len(value) == 1:
            datatype = f'single_{datatype}'
        else:
            datatype = f'multi_{datatype}'
        return value, datatype
    
    
    def _formulate_equal_query(self, last_val, last_val_name, value):
        qstr = "FILTER (?%s=\"%s\"^^xsd:%s)"%(last_val_name, 
                                              str(value[0]), 
                                              last_val['dtype'])
        return qstr
    
    def _formulate_or_string_query(self, last_val, last_val_name, value):
        qstr = "FILTER (?%s=\"%s\"^^xsd:%s || ?%s=\"%s\"^^xsd:%s)"%(last_val_name, 
                                                                    str(value[0]), 
                                                                    last_val['dtype'],
                                                                    last_val_name, 
                                                                    str(value[1]), 
                                                                    last_val['dtype'],)
        return qstr
    
    def _formulate_range_number_query(self, last_val, last_val_name, value):
        value = np.sort(value)
        qstr = "FILTER (?%s >= \"%s\"^^xsd:%s && ?%s <= \"%s\"^^xsd:%s)"%(last_val_name, 
                                                                    str(value[0]), 
                                                                    last_val['dtype'],
                                                                    last_val_name, 
                                                                    str(value[1]), 
                                                                    last_val['dtype'],)
        return qstr

        
            
            
            

        
        