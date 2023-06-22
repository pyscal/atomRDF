import numpy as np
import spglib

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

def get_lattice_vector(sys, cartesian=False):
    box = [[sys._structure_dict['scaling_factors'][0], 0, 0], 
           [0, sys._structure_dict['scaling_factors'][1], 0], 
           [0, 0, sys._structure_dict['scaling_factors'][2]]]
    if cartesian:
        box[0][0] = box[0][0]*sys.atoms._lattice_constant 
        box[1][1] = box[1][1]*sys.atoms._lattice_constant 
        box[2][2] = box[2][2]*sys.atoms._lattice_constant
    return box

def get_bravais_lattice(sys):
    lattice = sys.atoms._lattice 
    if lattice == "l12":
        lattice = "https://www.wikidata.org/wiki/Q3006714"
    elif lattice == "b2":
        lattice = "https://www.wikidata.org/wiki/Q851536"
    elif lattice == "diamond":
        lattice = "https://www.wikidata.org/wiki/Q3006714"
    elif lattice == "hcp":
        lattice = "https://www.wikidata.org/wiki/Q663314"
    elif lattice == "a15":
        lattice = "a15"
    elif lattice == "bcc":
        lattice = "https://www.wikidata.org/wiki/Q851536"
    elif lattice == "fcc":
        lattice = "https://www.wikidata.org/wiki/Q3006714"
    return lattice
    
def get_space_group(sys):
    box = get_lattice_vector(sys)
    direct_coordinates = sys._structure_dict['positions']
    atom_types = sys._structure_dict['species']
    results = spglib.get_symmetry_dataset((box,
    direct_coordinates, atom_types))
    return results["international"], results["number"]
    
def get_basis(sys):
    occ_numbers = sys._structure_dict['species']
    tdict = sys.atoms._type_dict
    vals = [val for key, val in tdict.items()]
    
    if vals[0] is not None:
        occ_numbers = [tdict[x] for x in occ_numbers]
    return occ_numbers

    
    

