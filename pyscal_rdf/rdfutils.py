import pyscal_rdf.properties as prp
import numpy as np

def convert_to_dict(sys):
    """
    Convert a pyscal System object to data dictionary
    
    Parameters
    ----------
    sys: pyscal System
        input system
    
    Returns
    -------
    info: dict
        dict with parsed information
    """
    info = {}
    
    
    #Simulation cell properties
    #--------------------------------------------------------------    
    info["ChemicalCompositionElement"] = list(sys.composition.keys())
    info["ChemicalCompositionRatio"] = [val for key, val in sys.composition.items()]    
    info["CellVolume"] = sys.volume
    info["NumberOfAtoms"] = sys.natoms
    
    box_dims = sys.box_dimensions
    info["SimulationCellLengthX"] = box_dims[0]
    info["SimulationCellLengthY"] = box_dims[1]
    info["SimulationCellLengthZ"] = box_dims[2]

    info["SimulationCellVectorA"] = sys.box[0]
    info["SimulationCellVectorB"] = sys.box[1]
    info["SimulationCellVectorC"] = sys.box[2]
    
    info["SimulationCellAngleAlpha"] = prp.get_angle(sys.box[0], sys.box[1])
    info["SimulationCellAngleBeta"] = prp.get_angle(sys.box[1], sys.box[2])
    info["SimulationCellAngleGamma"] = prp.get_angle(sys.box[2], sys.box[0])    
    
    #Atom properties
    #--------------------------------------------------
    if sys.atoms.species[0] is not None:
        info["Element"] = sys.atoms.species
    else:
        info["Element"] = sys.atoms.types
    info["Coordination"] = prp.get_coordination(sys)
    info["Positions"] = sys.atoms.positions
    
    #unit cell things
    #--------------------------------------------------
    #Unit cell properties
    info["LatticeParameter"] = sys.atoms._lattice_constant
    
    #do only for unit cell
    if sys._structure_dict is None:
        pass
        #maybe try to find space group with warnings!
    else:
        symbol, number = prp.get_space_group(sys)
        info["SpaceGroupSymbol"] = symbol
        info["SpaceGroupNumber"] = number
        info["CrystalStructureName"] = sys.atoms._lattice 
        info["BravaisLattice"] = prp.get_bravais_lattice(sys)
        info["BasisPositions"] = sys._structure_dict['positions']
        info["BasisOccupancy"] = prp.get_basis(sys)
        info["LatticeVectors"] = prp.get_lattice_vector(sys)
        
    return info