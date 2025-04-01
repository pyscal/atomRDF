import pyscal3.structure_creator as pcs
from atomrdf.system import System
from atomrdf.build.buildutils import _declass
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict

def bulk(
    element,
    structure=None,
    lattice_constant=1.00,
    repetitions=None,
    ca_ratio=1.633,
    noise=0,
    primitive=False,
    graph=None,
    names=False,
    label=None,
):
    """
    Create a crystal structure using the specified parameters.

    Parameters:
    -----------
    element : str or None, optional
        The element to use for the crystal. Default is None.
    structure : str
        The crystal structure to create.
    lattice_constant : float, optional
        The lattice constant of the crystal. Default is 1.00.
    repetitions : tuple or None, optional
        The number of repetitions of the crystal structure in each direction. Default is None.
    ca_ratio : float, optional
        The c/a ratio of the crystal. Default is 1.633.
    noise : float, optional
        The amount of noise to add to each atom position. Default is 0.
    primitive : bool, optional
        Whether to create a primitive cell. Default is False.
    graph : atomrdf.KnowledgeGraph, optional
        The graph object to use for the crystal. Default is None.
        The structure is added to the KnowledgeGraph only if this option is provided.
    names : bool, optional
        If provided, human-readable names will be assigned to each property. If False, random IDs will be used. Default is False.

    Returns:
    --------
    s : object
        The atomrdf.Structure object representing the generated crystal structure.
    """
    if repetitions is None:
        repetitions = [1, 1, 1]
    atoms, box, sdict = pcs.make_crystal(
        structure,
        lattice_constant = _declass(lattice_constant),
        repetitions=repetitions,
        ca_ratio = _declass(ca_ratio),
        noise=noise,
        element=element,
        return_structure_dict=True,
        primitive=primitive,
    )
    if 'repetitions' not in sdict.keys():
        sdict['repetitions'] = repetitions

    s = System(graph=graph, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = _declass(lattice_constant)
    s._structure_dict = sdict
    s.label = label
    s.to_graph()
    s.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    s.add_property_mappings(ca_ratio, mapping_quantity='lattice_constant')    
    return s

def lattice(
    positions,
    types,
    box,
    lattice_constant=1.00,
    repetitions=None,
    noise=0,
    element=None,
    graph=None,
    names=False,
    label=None,
):
    """
    Generate a general lattice structure.

    Parameters:
    -----------
    positions : array_like
        The atomic positions in the lattice.
    types : array_like
        The atomic types corresponding to the positions.
    box : array_like
        The box dimensions of the lattice.
    lattice_constant : float, optional
        The lattice constant, defaults to 1.00.
    repetitions : array_like, optional
        The number of repetitions of the lattice in each direction.
    noise : float, optional
        The amount of noise to add to the lattice positions, defaults to 0.
    element : str, optional
        The chemical elements associated with the atoms. Should be equal to the number of unique types.
    graph : atomrdf.KnowledgeGraph, optional
        The graph object to store the lattice structure, defaults to None.
        The structure is added to the KnowledgeGraph only if this option is provided.
    names : bool, optional
        If True, human readable names instead of random ids will be created. Default is False.

    Returns:
    --------
    s : object
        The atomrdf.Structure object representing the generated lattice structure.

    """
    atoms, box, sdict = pcs.general_lattice(
        positions,
        types,
        box,
        lattice_constant=_declass(lattice_constant),
        repetitions=repetitions,
        noise=noise,
        element=element,
        return_structure_dict=True,
    )

    if 'repetitions' not in sdict.keys():
        sdict['repetitions'] = repetitions

    s = System(graph=graph, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = "custom"
    s.atoms._lattice_constant = _declass(lattice_constant)
    s._structure_dict = sdict
    s.label = label
    s.to_graph()
    s.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    
    return s
ls
