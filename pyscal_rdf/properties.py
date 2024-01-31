import numpy as np
import spglib

# DATADICT properties
#------------------------------------------
bravais_lattice_dict = {
    "l12": "https://www.wikidata.org/wiki/Q3006714",
    "b2": "https://www.wikidata.org/wiki/Q851536",
    "diamond": "https://www.wikidata.org/wiki/Q3006714",
    "hcp": "https://www.wikidata.org/wiki/Q663314",
    "a15": "a15",
    "bcc": "https://www.wikidata.org/wiki/Q851536",
    "fcc": "https://www.wikidata.org/wiki/Q3006714",
}

# SIMCELL properties
#--------------------------------------------
def get_chemical_composition(system):
    return system.composition

def get_cell_volume(system):
    return system.volume

def get_number_of_atoms(system):
    return system.natoms

def get_simulation_cell_length(system):
    return system.box_dimensions

def get_simulation_cell_vector(system):
    return system.box

def get_simulation_cell_angle(system):
    return [_get_angle(system.box[0], system.box[1]),
            _get_angle(system.box[1], system.box[2]),
            _get_angle(system.box[2], system.box[0])]

# LATTICE properties
#--------------------------------------------

def get_lattice_angle(system):
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return [_get_angle(system._structure_dict["box"][0], system._structure_dict["box"][1]),
            _get_angle(system._structure_dict["box"][1], system._structure_dict["box"][2]),
            _get_angle(system._structure_dict["box"][2], system._structure_dict["box"][0])]
    else:
        return [None, None, None]

def get_lattice_parameter(system):
    if system.atoms._lattice_constant is None:
        return [None, None, None]
    else:
        if system._structure_dict is not None:
            if "box" in system._structure_dict.keys():
                return [np.linalg.norm(system._structure_dict["box"][0])*system.atoms._lattice_constant,
                        np.linalg.norm(system._structure_dict["box"][1])*system.atoms._lattice_constant,
                        np.linalg.norm(system._structure_dict["box"][2])*system.atoms._lattice_constant]
        return [system.atoms._lattice_constant, 
                system.atoms._lattice_constant, 
                system.atoms._lattice_constant]

def get_crystal_structure_name(system):
    if system._structure_dict is None:
        return None
    return system.atoms._lattice

def get_bravais_lattice(system):
    if system._structure_dict is None:
        return None
    if system.atoms._lattice in bravais_lattice_dict.keys():
        return bravais_lattice_dict[system.atoms._lattice]
    return None

def get_basis_positions(system):
    if system._structure_dict is None:
        return None
    if "positions" in system._structure_dict.keys():
        return system._structure_dict["positions"]
    return None

def get_basis_occupancy(system):
    if system._structure_dict is None:
        return None

    if "species" in system._structure_dict.keys():
        occ_numbers = system._structure_dict['species']
        tdict = system.atoms._type_dict
        vals = [val for key, val in tdict.items()]
        
        if vals[0] is not None:
            occ_numbers = [tdict[x] for x in occ_numbers]
        return occ_numbers
    return None

def get_lattice_vector(system):
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return system._structure_dict["box"]
    return [None, None, None]

def get_spacegroup_symbol(system):
    if system._structure_dict is None:
        return None
    try:
        results = _get_symmetry_dict(system)
        return results[0]
    except:
        return None

def get_spacegroup_number(system):
    if system._structure_dict is None:
        return None
    try:
        results = _get_symmetry_dict(system)
        return results[1]
    except:
        return None

# ATOM attributes
#--------------------------------------------
def get_position(system):
    return system.atoms.positions

def get_species(system):
    return system.atoms.species



# SUPPORT functions
#--------------------------------------------
def _get_angle(vec1, vec2):
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

def _get_symmetry_dict(system):
    box = get_lattice_vector(system)
    direct_coordinates = get_basis_positions(system)
    atom_types = system._structure_dict['species']

    results = spglib.get_symmetry_dataset((box,
    direct_coordinates, atom_types))
    return results["international"], results["number"]    
