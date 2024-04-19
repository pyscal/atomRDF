import numpy as np
import spglib

# DATADICT properties
# ------------------------------------------
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
# --------------------------------------------
def get_chemical_composition(system):
    """
    Get the chemical composition of the system.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    composition : dict
        A dictionary containing the chemical elements as keys and their corresponding counts as values.

    """
    return system.composition


def get_cell_volume(system):
    """
    Get the volume of the simulation cell.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    volume : float
        The volume of the simulation cell.

    """
    return system.volume


def get_number_of_atoms(system):
    """
    Get the number of atoms in the system.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    natoms : int
        The number of atoms in the system.

    """
    return system.natoms


def get_simulation_cell_length(system):
    """
    Get the length of the simulation cell.

    Parameters
    ----------
    system : object
        The system object.

    Returns
    -------
    length : list
        A list containing the length of each dimension of the simulation cell.

    """
    return system.box_dimensions


def get_simulation_cell_vector(system):
    """
    Get the simulation cell vector of the given system.

    Parameters
    ----------
    system : object
        The system object containing the simulation cell information.

    Returns
    -------
    numpy.ndarray
        The simulation cell vector of the system.

    """
    return system.box


def get_simulation_cell_angle(system):
    """
    Get the angles between the vectors of the simulation cell.

    Parameters
    ----------
    system : object
        The system object containing the simulation cell information.

    Returns
    -------
    angles : list
        A list containing the angles between the vectors of the simulation cell.

    """
    return [
        _get_angle(system.box[0], system.box[1]),
        _get_angle(system.box[1], system.box[2]),
        _get_angle(system.box[2], system.box[0]),
    ]


# LATTICE properties
# --------------------------------------------


def get_lattice_angle(system):
    """
    Calculate the lattice angles of a given system.

    Parameters
    ----------
    system : object
        The system object containing the structure information.

    Returns
    -------
    list
        A list of three lattice angles in degrees. If the structure information is not available, [None, None, None] is returned.

    """
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return [
            _get_angle(
                system._structure_dict["box"][0], system._structure_dict["box"][1]
            ),
            _get_angle(
                system._structure_dict["box"][1], system._structure_dict["box"][2]
            ),
            _get_angle(
                system._structure_dict["box"][2], system._structure_dict["box"][0]
            ),
        ]
    else:
        return [None, None, None]


def get_lattice_parameter(system):
    """
    Calculate the lattice parameters of a system.

    Parameters
    ----------
    system : object
        The system object containing information about the atoms and structure.

    Returns
    -------
    list
        A list containing the lattice parameters of the system. If the lattice constant is not available,
        [None, None, None] is returned. If the system structure is available, the lattice parameters are
        calculated based on the box dimensions. Otherwise, the lattice constant is returned for all three
        dimensions.

    Examples
    --------
    >>> system = System()
    >>> system.atoms._lattice_constant = 3.5
    >>> system._structure_dict = {"box": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}
    >>> get_lattice_parameter(system)
    [3.5, 3.5, 3.5]

    >>> system.atoms._lattice_constant = None
    >>> get_lattice_parameter(system)
    [None, None, None]
    """

    if system.atoms._lattice_constant is None:
        return [None, None, None]
    else:
        if system._structure_dict is not None:
            if "box" in system._structure_dict.keys():
                return [
                    np.linalg.norm(system._structure_dict["box"][0])
                    * system.atoms._lattice_constant,
                    np.linalg.norm(system._structure_dict["box"][1])
                    * system.atoms._lattice_constant,
                    np.linalg.norm(system._structure_dict["box"][2])
                    * system.atoms._lattice_constant,
                ]
        return [
            system.atoms._lattice_constant,
            system.atoms._lattice_constant,
            system.atoms._lattice_constant,
        ]


def get_crystal_structure_name(system):
    """
    Get the name of the crystal structure for a given system.

    Parameters
    ----------
    system : object
        The system object containing the crystal structure information.

    Returns
    -------
    str or None
        The name of the crystal structure if available, otherwise None.

    """
    if system._structure_dict is None:
        return None
    return system.atoms._lattice


def get_bravais_lattice(system):
    """
    Get the Bravais lattice of a given system.

    Parameters
    ----------
    system : object
        The system object for which the Bravais lattice is to be determined.

    Returns
    -------
    str or None
        The Bravais lattice of the system, or None if the system's structure dictionary is not available or the lattice is not found in the dictionary.

    """
    if system._structure_dict is None:
        return None
    if system.atoms._lattice in bravais_lattice_dict.keys():
        return bravais_lattice_dict[system.atoms._lattice]
    return None


def get_basis_positions(system):
    """
    Get the basis positions from the given system.

    Parameters
    ----------
    system : object
        The system object containing the structure dictionary.

    Returns
    -------
    numpy.ndarray or None
        The basis positions if available, otherwise None.
    """
    if system._structure_dict is None:
        return None
    if "positions" in system._structure_dict.keys():
        return system._structure_dict["positions"]
    return None


# def get_basis_occupancy(system):
#    if system._structure_dict is None:
#        return None

#    if "species" in system._structure_dict.keys():
#        occ_numbers = system._structure_dict['species']
#        tdict = system.atoms._type_dict
#        vals = [val for key, val in tdict.items()]

#        if vals[0] is not None:
#            occ_numbers = [tdict[x] for x in occ_numbers]
#        return occ_numbers
#    return None


def get_lattice_vector(system):
    """
    Get the lattice vector of a system.

    Parameters
    ----------
    system : object
        The system object containing the structure information.

    Returns
    -------
    list
        A list representing the lattice vector of the system. If the structure
        dictionary is not available or the lattice vector is not defined, it
        returns [None, None, None].
    """
    if system._structure_dict is None:
        return [None, None, None]
    if "box" in system._structure_dict.keys():
        return system._structure_dict["box"]
    return [None, None, None]


def get_spacegroup_symbol(system):
    """
    Get the symbol of the spacegroup for a given system.

    Parameters:
        system (object): The system object for which to retrieve the spacegroup symbol.

    Returns:
        str: The symbol of the spacegroup if available, otherwise None.
    """
    if system._structure_dict is None:
        return None
    try:
        results = _get_symmetry_dict(system)
        return results[0]
    except:
        return None


def get_spacegroup_number(system):
    """
    Get the spacegroup number of a given system.

    Parameters
    ----------
    system : object
        The system object for which the spacegroup number is to be determined.

    Returns
    -------
    int or None
        The spacegroup number of the system if it is available, otherwise None.
    """
    if system._structure_dict is None:
        return None
    try:
        results = _get_symmetry_dict(system)
        return results[1]
    except:
        return None


# ATOM attributes
# --------------------------------------------
def get_position(system):
    """
    Get the positions of the atoms in the system.

    Parameters
    ----------
    system : object
        The system object containing the atom positions.

    Returns
    -------
    numpy.ndarray or None
        The positions of the atoms if available, otherwise None.

    """
    return system.atoms.positions


def get_species(system):
    """
    Get the species of atoms in the given system.

    Parameters
    ----------
    system : System
        The system object containing atoms.

    Returns
    -------
    list
        A list of species of atoms in the system.

    """
    return system.atoms.species


# SUPPORT functions
# --------------------------------------------
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
    return np.round(
        np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
        * 180
        / np.pi,
        decimals=2,
    )


def _get_symmetry_dict(system):
    box = get_lattice_vector(system)
    direct_coordinates = get_basis_positions(system)
    atom_types = system._structure_dict["species"]

    results = spglib.get_symmetry_dataset((box, direct_coordinates, atom_types))
    return results["international"], results["number"]
