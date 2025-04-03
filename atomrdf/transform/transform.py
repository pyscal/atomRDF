import pyscal3.operations.operations as operations
from pyscal3.atoms import AttrSetter, Atoms
from atomrdf.structure import System
import numpy as np

def repeat(system, repetitions):
    """
    Repeat the system in each direction by the specified number of times.

    Parameters
    ----------
    repetitions : tuple
        The number of times to repeat the system in each direction.

    Returns
    -------
    None

    Notes
    -----
    The system is repeated in each direction by the specified number of times.
    """
    new_system = system.duplicate()
    new_system = operations.repeat(new_system, repetitions)
    if new_system._structure_dict is None:
        new_system._structure_dict = {}
    new_system._structure_dict["repetitions"] = repetitions
    new_system.to_graph()
    new_system.copy_defects(system.sample)
    return new_system

def rotate(sys, rotation_vectors, graph=None, label=None):
    try:
        from atomman.defect.Dislocation import Dislocation
        import atomman as am
        import atomman.unitconvert as uc
    except ImportError:
        raise ImportError("This function requires the atomman package to be installed")

    box = am.Box(
        avect=sys.box[0],
        bvect=sys.box[1],
        cvect=sys.box[2],
    )
    
    atoms = am.Atoms(
        atype=sys.atoms.types, pos=sys.atoms.positions
    )
    
    element = [val for key, val in sys.atoms._type_dict.items()]

    system = am.System(
        atoms=atoms, 
        box=box, 
        pbc=[True, True, True], 
        symbols=element, 
        scale=False,
    )

    #now rotate with atomman
    system = system.rotate(rotation_vectors)

    #now convert back and return the system
    box = [system.box.avect, 
        system.box.bvect, 
        system.box.cvect]
    
    atom_df = system.atoms_df()
    types = [int(x) for x in atom_df.atype.values]
    
    species = []
    for t in types:
        species.append(element[int(t) - 1])

    positions = np.column_stack(
        (atom_df["pos[0]"].values, 
        atom_df["pos[1]"].values, 
        atom_df["pos[2]"].values)
    )

    atom_dict = {"positions": positions, "types": types, "species": species}
    atom_obj = Atoms()
    atom_obj.from_dict(atom_dict)
    
    output_structure = System()
    output_structure.box = box
    output_structure.atoms = atom_obj
    #output_structure = output_structure.modify.remap_to_box()
    if graph is not None:
        output_structure.graph = graph
    else:
        output_structure.graph = sys.graph
    output_structure.atoms._lattice = sys.atoms._lattice
    output_structure.atoms._lattice_constant = sys.atoms._lattice_constant
    output_structure._structure_dict = sys._structure_dict
    if label is not None:
        output_structure.label = label
    else:
        output_structure.label = sys.label
    output_structure.to_graph()
    output_structure.copy_defects(sys.sample)
    if output_structure.graph is not None:
        sys.add_rotation_triples(rotation_vectors, output_structure.sample)
    return output_structure

def translate(system,
                translation_vector, 
                plane=None, distance=None, 
                reverse_orientation=False, 
                copy_structure=True,
                add_triples=True):
    original_sample = system.sample
    if copy_structure:
        sys = system.duplicate()
        #and add this new structure to the graph
        sys.to_graph()
        sys.copy_defects(system.sample)
    else:
        sys = system

    if plane is not None:
        if distance is None:
            raise ValueError('distance needs to be provided')
    
    if plane is not None:
        sys.select_by_plane(plane, distance, reverse_orientation=reverse_orientation)

    if not len(translation_vector) == 3:
        raise ValueError("translation vector must be of length 3")
    
    translation_vector = np.array(translation_vector)

    for x in range(len(sys.atoms['positions'])):
        if sys.atoms['condition'][x]:
            sys.atoms['positions'][x] += translation_vector
    
    if plane is not None:
        sys.remove_selection()

    if add_triples:
        sys.add_translation_triples(translation_vector, plane, distance, original_sample)
    return sys
