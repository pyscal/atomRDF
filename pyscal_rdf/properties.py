import numpy as np
from pyscal.operations.symmetry import Symmetry

def get_angle(vec1, vec2):
    """
    Get angle between two vectors in degrees
    
    Parameters
    ----------
    vec1: list
        first vector
    
    vec2: list
        second vector
    
    Returns
    -------
    angle: float
        angle in degrees
    
    Notes
    -----
    Angle is rounded to two decimal points
    
    """
    return np.round(np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))*180/np.pi, decimals=2)

def get_coordination(sys):
    sys.find_neighbors(method="cutoff")
    coordination = [len(x) for x in sys.atoms.neighbors.index]
    return coordination


def get_space_group(sys):
    sym = Symmetry(sys)
    sym.calculate()
    return sym.output.international_space_group_number, sym.output.international_symbol
    

