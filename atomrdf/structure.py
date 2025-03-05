"""
This module provides functions for creating and manipulating atomic structures. It includes functions for creating crystals,
general lattices, dislocations, grain boundaries, and reading structures from files. The module also provides functionality for
adding interstitial impurities, substituting atoms, deleting atoms, and adding vacancies.
The structures can be converted to RDF graphs using the atomrdf library.
The main object in this module is the System class, which extends the functionality of the pyscal3.core.System class and provides additional methods for working with atomic structures.
"""

import numpy as np
import copy
from functools import partial, update_wrapper
import os
import yaml
import uuid
import json
import shutil
import tarfile
import warnings
from ase.io import write
from ase.build import cut

import pyscal3.structure_creator as pcs
from pyscal3.grain_boundary import GrainBoundary
from pyscal3.atoms import AttrSetter, Atoms
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict
from pyscal3.formats.ase import convert_snap

import pyscal3.operations.input as inputmethods
import pyscal3.operations.serialize as serialize
import pyscal3.operations.visualize as visualize
import pyscal3.operations.operations as operations

import atomrdf.json_io as json_io
import atomrdf.properties as prp
from atomrdf.sample import Property
import atomrdf.io as aio

from rdflib import Graph, Namespace, XSD, RDF, RDFS, BNode, URIRef
from atomrdf.namespace import CMSO, LDO, PLDO, PODO, UNSAFEASMO, UNSAFECMSO, PROV, Literal, ASMO

# read element data file
file_location = os.path.dirname(__file__).split("/")
file_location = "/".join(file_location[:-1])
file_location = os.path.join(os.path.dirname(__file__), "data/element.yml")
with open(file_location, "r") as fin:
    element_indetifiers = yaml.safe_load(fin)

#declassing special variables
def _declass(item):
    if isinstance(item, Property):
        return item.value
    else:
        return item
        
def _make_crystal(
    structure,
    lattice_constant=1.00,
    repetitions=None,
    ca_ratio=1.633,
    noise=0,
    element=None,
    primitive=False,
    graph=None,
    names=False,
    label=None,
):
    """
    Create a crystal structure using the specified parameters.

    Parameters:
    -----------
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
    element : str or None, optional
        The element to use for the crystal. Default is None.
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


def _make_general_lattice(
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

def _make_stacking_fault(
    slip_plane,
    displacement_a,
    displacement_b=0,
    slip_direction_a=None,
    slip_direction_b=None,
    vacuum=0,
    minwidth=15,
    even=True,
    minimum_r=None,
    relative_fault_position=0.5,
    structure=None,
    element=None,
    lattice_constant=1.00,
    repetitions=None,
    ca_ratio=1.633,
    noise=0,
    primitive=False,
    graph=None,
    names=False,
    label=None,
    return_atomman_dislocation=False,
):
    """
    Generate a stacking fault structure.

    Parameters
    ----------
    slip_system : list of lists, shape (2 x 3) or (2 x 4)
        the slip system for the given system. The input should of type [[u, v, w], [h, k, l]].
        [u, v, w] is the slip direction and [h, k, l] is the slip plane.

        For HCP systems, the input should be [[u, v, w, z], [h, k, l, m]].
    
    distance : float
        Distance for translating one half of the cell along the [h k l] direction. Default is 1.
    """
    try:
        import atomman as am
        import atomman.unitconvert as uc
    except ImportError:
        raise ImportError("This function requires the atomman package to be installed")

    
    #we are good to go now
    #we proceed with the atomman code
    if element is None:
        raise ValueError("Please provide element")
    if element in element_dict.keys():
        structure = element_dict[element]["structure"]
        lattice_constant = element_dict[element]["lattice_constant"]
    else:
        raise ValueError("Please provide structure")
    if structure == "hcp":
        if len(slip_plane) != 4:
            raise ValueError("For hcp systems, slip plane should be of length 4")
        #routine for hcp
        a = uc.set_in_units(_declass(lattice_constant), 'angstrom')
        c = uc.set_in_units(_declass(ca_ratio), 'angstrom')
        atoms = am.Atoms(pos=[[0.0, 0.0, 0.0], [1/3, 2/3, 0.5]])        
        ucell = am.System(atoms=atoms, box=am.Box.hexagonal(a, c), scale=True, symbols=element)
        sd = structure_dict["hcp"]["primitive"]
    else:
        #routine for others
        #extract the structure vectors
        if primitive:
            sd = structure_dict[structure]["primitive"]
        else:
            sd = structure_dict[structure]["conventional"]
        vectors = _declass(lattice_constant)*np.array(sd["box"])
        #positions
        positions = np.array(sd["positions"])
        types = np.array(sd["species"])
        
        # create a structure with the info
        box = am.Box(
            avect=vectors[0],
            bvect=vectors[1],
            cvect=vectors[2],
        )

        atoms = am.Atoms(
            atype=types, pos=positions,
        )
        ucell = am.System(atoms=atoms, box=box, scale=True, symbols=element)
    
    #seupersize
    if repetitions is not None:
        if isinstance(repetitions, int):
            repetitions = [repetitions, repetitions, repetitions]
        ucell = ucell.supersize(*repetitions)

    sf = am.defect.StackingFault(slip_plane, ucell)
    if slip_direction_a is not None:
        sf.a1vect_uvw = slip_direction_a
    if slip_direction_b is not None:
        sf.a2vect_uvw = slip_direction_b
    surfacesystem = sf.surface(shift=sf.shifts[0], 
                            minwidth=minwidth, 
                            even=even,
                            vacuumwidth=vacuum)
    if relative_fault_position != 0.5:
        sf.faultpos_rel = relative_fault_position
    faultsystem = sf.fault(a1=displacement_a, 
                        a2=displacement_b)

    #get displacements
    displ = am.displacement(surfacesystem, faultsystem)

    box = [faultsystem.box.avect, faultsystem.box.bvect, faultsystem.box.cvect]
    atom_df = faultsystem.atoms_df()
    types = [int(x) for x in atom_df.atype.values]
    if element is not None:
        if isinstance(element, str):
            element = [element]
        species = []
        for t in types:
            species.append(element[int(t) - 1])
    else:
        species = [None for x in range(len(types))]

    positions = np.column_stack(
        (atom_df["pos[0]"].values, atom_df["pos[1]"].values, atom_df["pos[2]"].values)
    )

    positions = positions + displ

    #create datadict
    datadict = {}
    datadict["lattice"] = structure
    datadict["lattice_constant"] = _declass(lattice_constant)
    datadict["box"] = sd["box"]
    datadict["positions"] = sd["positions"]
    datadict["repetitions"] = repetitions

    output_structure = System()
    output_structure.box = box
    atoms = Atoms()
    atoms.from_dict({"positions": positions, "species": species, "types": types})
    output_structure.atoms = atoms
    #output_structure = output_structure.modify.remap_to_box()
    output_structure.lattice_properties = datadict
    output_structure.label = label
    output_structure.graph = graph
    output_structure.to_graph()
    output_structure.add_stacking_fault({"plane":slip_plane, "displacement":sf.a1vect_uvw})
    output_structure.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    output_structure.add_property_mappings(ca_ratio, mapping_quantity='lattice_constant')

    if return_atomman_dislocation:
        return output_structure, sf, surfacesystem, faultsystem
    return output_structure


def _make_dislocation(
    slip_system,
    dislocation_line,
    elastic_constant_dict,
    burgers_vector=None,
    dislocation_type="monopole",
    structure=None,
    element=None,
    lattice_constant=1.00,
    repetitions=None,
    ca_ratio=1.633,
    noise=0,
    primitive=False,
    graph=None,
    names=False,
    label=None,
    return_atomman_dislocation=False,
):
    """
    Generate a dislocation structure. Wraps the atomman.defect.Dislocation class.

    Parameters
    ----------
    slip_system : list of lists, shape (2 x 3)
        the slip system for the given system. The input should of type [[u, v, w], [h, k, l]].
        [u, v, w] is the slip direction and [h, k, l] is the slip plane.
    dislocation_line : numpy array of length 3
        The dislocation line direction.
        This determines the type of dislocation to generate, screw, edge or mixed.
    burgers_vector : scalar or numpy array of length 3, optional
        if a scalar value (b) is provided, the burgers vector is assumed to be b*[u, v, w].
        if a numpy array is provided, the burgers vector is set to this value.
        Default is equal to slip direction [u, v, w] from slip_system.
    elastic_constant_dict : dict
        Dictionary of elastic constants. The keys should be in Voigt notation.
    burgers_vector : numpy array of length 3
        The Burgers vector of the dislocation.
    dislocation_type : str, optional
        The type of dislocation to generate. Default is "monopole".
    structure : crystal lattice to be used
        The crystal structure to use to create the bulk system. Either structure or element must be given.
    element : str, optional
        The atomic symbol according to which the bulk structure is to be created. Either structure or element must be given.
    lattice_constant : float, optional
        The lattice constant to use for the generated crystal structure. Default is 1.00.
    repetitions : tuple, optional
        The number of times to repeat the unit cell in each of the three Cartesian directions. Default is None.
    ca_ratio : float, optional
        The c/a ratio to use for the generated crystal structure. Default is 1.633. Used only if structure is hcp.
    noise : float, optional
        Magnitude of random noise to add to the atomic positions. Default is 0.
    primitive : bool, optional
        If True, the generated crystal structure will be converted to a primitive cell. Default is False.
    graph :atomrdf.KnowledgeGraph, optional
        A graph object representing the crystal structure. Default is None.
    names : bool, optional
        If True, the returned System will have atom names assigned based on the element and index. Default is False.

    Returns
    -------
    output_structure : atomrdf.System
        The generated dislocation structure.

    Raises
    ------
    ValueError
        If neither structure nor element is provided.

    Notes
    -----
    This function requires the atomman Python package to be installed.

    The elastic_constant_dict parameter should be a dictionary of elastic constants with keys corresponding to the
    following Voigt notation: "C11", "C12", "C13", "C14", "C15", "C16", "C22", "C23", "C24", "C25", "C26", "C33", "C34",
    "C35", "C36", "C44", "C45", "C46", "C55", "C56", "C66". The values should be given in GPa.

    The dislocation_type parameter can be set to "monopole" or "periodicarray". If set to "monopole", a single dislocation
    will be generated. If set to "periodicarray", a periodic array of dislocations will be generated.

    Needs atomman.
    """

    try:
        from atomman.defect.Dislocation import Dislocation
        import atomman as am
        import atomman.unitconvert as uc
    except ImportError:
        raise ImportError("This function requires the atomman package to be installed")
    
    slip_direction = slip_system[0]
    slip_plane = slip_system[1]
    if burgers_vector is None:
        burgers_vector = slip_direction
    elif np.isscalar(burgers_vector):
        burgers_vector = burgers_vector * np.array(slip_direction)
    elif len(burgers_vector) != 3:
        raise ValueError('burgers vector should be None, scalar, or of length 3')


    if structure is not None:
        # create a structure with the info
        input_structure = _make_crystal(
            structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=noise,
            element=element,
            primitive=primitive,
        )
    elif element is not None:
        if element in element_dict.keys():
            structure = element_dict[element]["structure"]
            lattice_constant = element_dict[element]["lattice_constant"]
        else:
            raise ValueError("Please provide structure")
        input_structure = _make_crystal(
            structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=noise,
            element=element,
            primitive=primitive,
        )
    else:
        raise ValueError("Provide either structure or element")

    for key, val in elastic_constant_dict.items():
        elastic_constant_dict[key] = uc.set_in_units(val, "GPa")
    C = am.ElasticConstants(**elastic_constant_dict)

    box = am.Box(
        avect=input_structure.box[0],
        bvect=input_structure.box[1],
        cvect=input_structure.box[2],
    )
    atoms = am.Atoms(
        atype=input_structure.atoms.types, pos=input_structure.atoms.positions
    )
    system = am.System(
        atoms=atoms, box=box, pbc=[True, True, True], symbols=element, scale=False
    )

    disc = Dislocation(
        system,
        C,
        burgers_vector,
        dislocation_line,
        slip_plane,
    )
    if dislocation_type == "monopole":
        disl_system = disc.monopole()
    elif dislocation_type == "periodicarray":
        disl_system = disc.periodicarray()

    box = [disl_system.box.avect, disl_system.box.bvect, disl_system.box.cvect]
    atom_df = disl_system.atoms_df()
    types = [int(x) for x in atom_df.atype.values]
    if element is not None:
        species = []
        for t in types:
            species.append(element[int(t) - 1])
    else:
        species = [None for x in range(len(types))]

    positions = np.column_stack(
        (atom_df["pos[0]"].values, atom_df["pos[1]"].values, atom_df["pos[2]"].values)
    )

    #find dislocation character
    angle = np.dot(dislocation_line, burgers_vector)/(np.linalg.norm(dislocation_line)*np.linalg.norm(burgers_vector))
    angle_rad = np.arccos(angle)
    angle_deg = np.degrees(angle_rad)

    disl_dict = {
        'BurgersVector': burgers_vector,
        'SlipPlane': slip_plane,
        'SlipDirection': slip_direction,
        'DislocationLine': dislocation_line,
        'DislocationCharacter': angle_deg,
    }

    # here we dont add repetitions, since we cannot guarantee
    atom_dict = {"positions": positions, "types": types, "species": species}
    atom_obj = Atoms()
    atom_obj.from_dict(atom_dict)
    output_structure = System()
    output_structure.box = box
    output_structure.atoms = atom_obj
    output_structure = output_structure.modify.remap_to_box()
    output_structure.label = label
    output_structure.graph = graph
    output_structure.to_graph()
    output_structure.add_dislocation(disl_dict)
    output_structure.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    output_structure.add_property_mappings(ca_ratio, mapping_quantity='lattice_constant')

    if return_atomman_dislocation:
        return output_structure, disc
    return output_structure

def _make_grain_boundary(
    axis,
    sigma,
    gb_plane,
    structure=None,
    element=None,
    lattice_constant=1,
    ca_ratio=1.633,
    repetitions=(1, 1, 1),
    overlap=0.0,
    gap=0.0,
    vacuum=0.0,
    delete_layer="0b0t0b0t",
    tolerance=  0.25,
    primitive=False,
    uc_a=1,
    uc_b=1,
    graph=None,
    names=False,
    label=None,
    backend='aimsgb',
    add_extras=False,         
):
    """
    Create a grain boundary system. GB can be created either with AIMSGB or GBCode.

    Parameters:
    -----------
    axis : tuple or list
        The rotation axis of the grain boundary.
        Used with backend 'aimsgb' and 'gbcode'.
    sigma : int
        The sigma value of the grain boundary.
        Used with backend 'aimsgb' and 'gbcode'.
    gb_plane : tuple or list
        The Miller indices of the grain boundary plane.
        Used with backend 'aimsgb' and 'gbcode'.
    backend : str, optional
        The backend to use to create the grain boundary. Default is 'aimsgb'.
        Some keyword arguments are only suitable for some backend.
    structure : the lattice structure to be used to create the GB, optional
        The lattice structure to populate the grain boundary with.
        Used with backend 'aimsgb' and 'gbcode'.
    element : str, optional
        The element symbol to populate the grain boundary with.
        Used with backend 'aimsgb' and 'gbcode'.
    lattice_constant : float, optional
        The lattice constant of the structure.
        Used with backend 'aimsgb' and 'gbcode'.
    repetitions : tuple or list, optional
        The number of repetitions of the structure that will be used to create the GB.
        Used only with 'gbcode'.
        For example, if (2,3,4) is provided, each grain will have these repetitions in (x,y,z) directions.
        For similar functionality in 'aimsgb', use 'uc_a' and 'uc_b'.
    overlap : float, optional
        The overlap between adjacent grain boundaries.
        Used only with 'gbcode'.
    vaccum : float, optional
        Adds space between the grains at one of the two interfaces
        that must exist due to periodic boundary conditions.
        Used only with 'aimsgb'.
    gap: float, optional
        Adds space between the grains at both of the two interfaces
        that must exist due to periodic boundary conditions.
        Used only with 'aimsgb'.
    delete_layer: str, optional
        To delete layers of the GB.
        Used only with 'aimsgb'.
    tolerance: float, optional
        Tolerance factor (in distance units) to determine whether two atoms
        are in the same plane.
        Used only with 'aimsgb'.
    primitive: bool, optional
        To generate primitive or non-primitive GB structure.
        Used only with 'aimsgb'.
    uc_a: int, optional
        Number of unit cells of left grain.
        Used only with 'aimsgb'.
    uc_b: int, optional
        Number of unit cells of right grain.
        Used only with 'aimsgb'.
    graph : atomrdf.KnowledgeGraph, optional
        The graph object to store the system.
        The system is only added to the KnowledgeGraph  if this option is provided.
    names : bool, optional
        If True human readable names will be assigned to each property. If False random ids will be used. Default is False.
    label: str, optional
        Add a label to the structure
    add_extras: bool, optional
        returns internal objects of the GB creation process.

    Returns:
    --------
    atomrdf.System
        The grain boundary system.

    Notes
    -----
    This function requires the aimsgb and pymatgen packages to be installed to use the 'aimsgb' backend.

    `repetitions` is used only with the 'gbcode' backend. 
    For similar functionality in 'aimsgb', use `uc_a` and `uc_b`. However, repetition in the third direction
    is not supported in 'aimsgb'. For a similar effect, after reaching the GB, `system.modify.repeat` function
    could be used with (1, 1, u_c).

    If 'gbcode' is used as backend, the specific type of GB is determined using the `find_gb_character` function
    When backend 'aimsgb' is used, this is attempted. If the type could not be found, a normal GB will be added in the annotation.

    """ 
    if backend == 'aimsgb':
        return _make_grain_boundary_aimsgb(
            axis,
            sigma,
            gb_plane,
            structure=structure,
            element=element,
            lattice_constant=lattice_constant,
            ca_ratio=ca_ratio,
            repetitions=repetitions,
            gap=gap,
            vacuum=vacuum,
            delete_layer=delete_layer,
            tolerance=  tolerance,
            primitive=primitive,
            uc_a=uc_a,
            uc_b=uc_b,
            graph=graph,
            names=names,
            label=label, 
            add_extras=add_extras,              
        )
    else:
        return _make_grain_boundary_gbcode(
            axis,
            sigma,
            gb_plane,
            structure=structure,
            element=element,
            lattice_constant=lattice_constant,
            repetitions=repetitions,
            overlap=overlap,
            graph=graph,
            names=names,
            label=label,
            add_extras=add_extras,
        )

def _make_grain_boundary_aimsgb(
    axis,
    sigma,
    gb_plane,
    structure=None,
    element=None,
    lattice_constant=1,
    ca_ratio=1.633,
    repetitions=(1, 1, 1),
    gap=0.0,
    vacuum=0.0,
    delete_layer="0b0t0b0t",
    tolerance=  0.25,
    primitive=False,
    uc_a=1,
    uc_b=1,
    graph=None,
    names=False,
    label=None,
    add_extras=False,  
): 
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
        from aimsgb import GrainBoundary as AIMSGrainBoundary
        from aimsgb import Grain as AIMSGrain
    except ImportError:
        raise ImportError("This function requires the aimsgb and pymatgen packages to be installed")
    

    if structure is not None:
        # create a structure with the info
        init_sys = _make_crystal(
            structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=0,
            element=element,
            primitive=primitive,
        )
    elif element is not None:
        if element in element_dict.keys():
            structure = element_dict[element]["structure"]
            lattice_constant = element_dict[element]["lattice_constant"]
        else:
            raise ValueError("Please provide structure")
        init_sys = _make_crystal(
            structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=0,
            element=element,
            primitive=primitive,
        )
    else:
        raise ValueError("Provide either structure or element")

    sdict = copy.deepcopy(init_sys._structure_dict)

    asesys = init_sys.write.ase()
    pmsys = AseAtomsAdaptor().get_structure(atoms=asesys)
    grain = AIMSGrain(pmsys.lattice, pmsys.species, pmsys.frac_coords)
    gb = AIMSGrainBoundary(axis=axis, sigma=sigma, 
                    plane=gb_plane, 
                    initial_struct=grain, 
                    uc_a=uc_a, 
                    uc_b=uc_b)
    gb_struct = AIMSGrain.stack_grains(
                grain_a = gb.grain_a,
                grain_b = gb.grain_b,
                vacuum = vacuum,
                gap=gap,
                direction = gb.direction,
                delete_layer=delete_layer,
                tol=tolerance,
                to_primitive=primitive,
            )
    asestruct = AseAtomsAdaptor().get_atoms(structure=gb_struct)
    sys = System.read.ase(asestruct, graph=None, names=names, label=label)
    sys.atoms._lattice = structure
    sys.atoms._lattice_constant = _declass(lattice_constant)
    sys._structure_dict = sdict
    sys.label = label
    sys.graph = graph
    sys.to_graph()
    sys.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    sys.add_property_mappings(ca_ratio, mapping_quantity='lattice_constant')

    try:
        gb_inb = GrainBoundary()
        gb_inb.create_grain_boundary(axis=axis, sigma=sigma, gb_plane=gb_plane)
        gb_type = gb_inb.find_gb_character()
    except:
        gb_type = None

    gb_dict = {
        "GBPlane": " ".join(np.array(gb_plane).astype(str)),
        "RotationAxis": axis,
        "MisorientationAngle": gb.theta[0],
        "GBType": gb_type,
        "sigma": gb.sigma,
    }
    sys.add_gb(gb_dict)

    if add_extras:
        sys._ase_system = asestruct
        sys._gb = gb
        sys._gb_struct = gb_struct
    return sys



def _make_grain_boundary_gbcode(
    axis,
    sigma,
    gb_plane,
    structure=None,
    element=None,
    lattice_constant=1,
    repetitions=(1, 1, 1),
    overlap=0.0,
    graph=None,
    names=False,
    label=None,
    add_extras=False,  
):
    """
    Create a grain boundary system.

    Parameters:
    -----------
    axis : tuple or list
        The rotation axis of the grain boundary.
    sigma : int
        The sigma value of the grain boundary.
    gb_plane : tuple or list
        The Miller indices of the grain boundary plane.
    structure : the lattice structure to be used to create the GB, optional
        The lattice structure to populate the grain boundary with.
    element : str, optional
        The element symbol to populate the grain boundary with.
    lattice_constant : float, optional
        The lattice constant of the structure.
    repetitions : tuple or list, optional
        The number of repetitions of the grain boundary structure in each direction.
    overlap : float, optional
        The overlap between adjacent grain boundaries.
    graph : atomrdf.KnowledgeGraph, optional
        The graph object to store the system.
        The system is only added to the KnowledgeGraph  if this option is provided.
    names : bool, optional
        If True human readable names will be assigned to each property. If False random ids will be used. Default is False.

    Returns:
    --------
    atomrdf.System
        The grain boundary system.

    """
    gb = GrainBoundary()
    gb.create_grain_boundary(axis=axis, sigma=sigma, gb_plane=gb_plane)

    if structure is not None:
        atoms, box, sdict = gb.populate_grain_boundary(
            structure,
            repetitions=repetitions,
            lattice_parameter=_declass(lattice_constant),
            overlap=overlap,
        )
    elif element is not None:
        atoms, box, sdict = gb.populate_grain_boundary(
            element, repetitions=repetitions, overlap=overlap
        )

    if 'repetitions' not in sdict.keys():
        sdict['repetitions'] = repetitions
    
    s = System(graph=None, names=names)
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = _declass(lattice_constant)
    s._structure_dict = sdict
    s.label = label
    s.graph = graph
    s.to_graph()
    s.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    
    gb_dict = {
        "GBPlane": " ".join(np.array(gb_plane).astype(str)),
        "RotationAxis": axis,
        "MisorientationAngle": gb.theta,
        "GBType": gb.find_gb_character(),
        "sigma": gb.sigma,
    }
    s.add_gb(gb_dict)
    return s


def _read_structure(
    filename,
    format="lammps-dump",
    graph=None,
    names=False,
    species=None,
    lattice=None,
    lattice_constant=None,
    basis_box=None,
    basis_positions=None,
    label=None,
    repetitions=None,
):
    """
    Read in structure from file or ase object

    Parameters
    ----------
    filename: string
        name of file

    format: string, optional
        format of the file

    graph: atomrdf.KnowledgeGraph, optional
        if provided, the structure will be added to the graph

    names: bool, optional
        if True, human readable names instead of random ids will be created.

    species: list, optional
        if provided lammps types will be matched to species. For example, if types 1 and 2 exist
        in the input file, and species = ['Li', 'Al'] is given, type 1 will be matched to 'Li' and
        type 2 will be matched to 'Al'

    lattice: str, optional
        currently supported lattices are {simple_cubic, bcc, fcc, hcp, dhcp, diamond, a15, l12, b2}
        if provided, information such as the unit cell, basis positions, space groups, etc are automatically added

    lattice_constant: float, optional
        specify the lattice constant of the system

    basis_box: 3x3 list, optional
        specify the basis unit cell. Not required if lattice is provided

    basis_positions: nX3 list, optional
        specify the relative positions of atoms in the unit cell. Not required if lattice is provided
    
    repetitions: tuple, optional
        specify the number of repetitions of the unit cell in each direction. Default is None.

    Returns
    -------
    atomrdf.System
    """
    datadict = {}
    if lattice is not None:
        if lattice in structure_dict.keys():
            datadict = structure_dict[lattice]["conventional"]
        datadict["lattice"] = lattice
    if lattice_constant is not None:
        datadict["lattice_constant"] = _declass(lattice_constant)
    if basis_box is not None:
        datadict["box"] = basis_box
    if basis_positions is not None:
        datadict["positions"] = basis_positions
    if repetitions is not None:
        datadict["repetitions"] = repetitions

    s = System(
        filename,
        format=format,
        species=species,
        graph=graph,
        names=names,
        warn_read_in=False,
    )
    s.lattice_properties = datadict
    s.label = label
    s.to_graph()
    s.add_property_mappings(lattice_constant, mapping_quantity='lattice_constant')
    return s


class System(pc.System):

    create = AttrSetter()
    # create.head = pcs
    mapdict = {}
    mapdict["lattice"] = {}
    for key in structure_dict.keys():
        mapdict["lattice"][key] = update_wrapper(
            partial(_make_crystal, key), _make_crystal
        )
    mapdict["lattice"]["custom"] = _make_general_lattice

    mapdict["element"] = {}
    for key in element_dict.keys():
        mapdict["element"][key] = update_wrapper(
            partial(
                _make_crystal,
                element_dict[key]["structure"],
                lattice_constant=element_dict[key]["lattice_constant"],
                element=key,
            ),
            pcs.make_crystal,
        )

    mapdict["defect"] = {}
    mapdict["defect"]["grain_boundary"] = _make_grain_boundary
    mapdict["defect"]["dislocation"] = _make_dislocation
    mapdict["defect"]["stacking_fault"] = _make_stacking_fault
    create._add_attribute(mapdict)

    read = AttrSetter()
    mapdict = {}
    mapdict["file"] = _read_structure
    mapdict["ase"] = update_wrapper(
        partial(_read_structure, format="ase"), _read_structure
    )
    read._add_attribute(mapdict)

    def __init__(
        self,
        filename=None,
        format="lammps-dump",
        compressed=False,
        customkeys=None,
        species=None,
        source=None,
        graph=None,
        names=False,
        warn_read_in=True,
    ):

        if (filename is not None) and warn_read_in:
            warnings.warn(
                "To provide additional information, use the System.read.file method"
            )

        super().__init__(
            filename=filename,
            format=format,
            compressed=compressed,
            customkeys=customkeys,
            species=species,
        )

        # this is the sample which will be stored
        self.sample = None
        self.label = None
        # the graph object should also be attached
        # for post-processing of structures
        self.graph = graph
        self.names = names
        self._material = None
        self._name = None
        self._atom_ids = None
        if source is not None:
            self.__dict__.update(source.__dict__)

        self.write = AttrSetter()
        mapdict = {}
        mapdict['ase'] = update_wrapper(partial(self.to_file, format='ase'), self.to_file)
        mapdict['pyiron'] = update_wrapper(partial(self.to_file, format='pyiron'), self.to_file)
        mapdict['quantum_espresso'] = update_wrapper(partial(self.to_file, format='quantum-espresso'), self.to_file)
        mapdict['file'] = self.to_file
        self.write._add_attribute(mapdict)

        self.show = AttrSetter()
        mapdict = {}
        mapdict['all'] = update_wrapper(partial(self._plot_system, plot_style='all'), self._plot_system)
        mapdict['continuous_property'] = update_wrapper(partial(self._plot_system, plot_style='continuous_property'), self._plot_system)
        mapdict['boolean_property'] = update_wrapper(partial(self._plot_system, plot_style='boolean_property'), self._plot_system)
        mapdict['selection'] = update_wrapper(partial(self._plot_system, plot_style='selection'), self._plot_system)
        self.show._add_attribute(mapdict)

        self.modify = AttrSetter()
        mapdict = {}
        mapdict["repeat"] = self.repeat
        mapdict["delete"] = self.delete
        mapdict["transform_to_cubic_cell"] = update_wrapper(partial(operations.extract_cubic_representation, self), operations.extract_cubic_representation)
        mapdict["remap_to_box"] = update_wrapper(partial(operations.remap_to_box, self), operations.remap_to_box)
        mapdict["remap_position_to_box"] = update_wrapper(partial(operations.remap_position_to_box, self), operations.remap_position_to_box)
        mapdict["embed_in_cubic_box"] = update_wrapper(partial(operations.embed_in_cubic_box, self), operations.embed_in_cubic_box)
        self.modify._add_attribute(mapdict)

        self.schema = AttrSetter()
        mapdict = {
            "material": {
                "element_ratio": partial(prp.get_chemical_composition, self),
                "crystal_structure": {
                    "name": partial(prp.get_crystal_structure_name, self),
                    "spacegroup_symbol": partial(prp.get_spacegroup_symbol, self),
                    "spacegroup_number": partial(prp.get_spacegroup_number, self),
                    "unit_cell": {
                        "bravais_lattice": partial(prp.get_bravais_lattice, self),
                        "lattice_parameter": partial(prp.get_lattice_parameter, self),
                        "angle": partial(prp.get_lattice_angle, self),
                    },
                },
            },
            "simulation_cell": {
                "volume": partial(prp.get_cell_volume, self),
                "number_of_atoms": partial(prp.get_number_of_atoms, self),
                "length": partial(prp.get_simulation_cell_length, self),
                "vector": partial(prp.get_simulation_cell_vector, self),
                "angle": partial(prp.get_simulation_cell_angle, self),
                "repetitions": partial(prp.get_repetitions, self),
            },
            "atom_attribute": {
                "position": partial(prp.get_position, self),
                "species": partial(prp.get_species, self),
            },
        }

        self.schema._add_attribute(mapdict)

    def _plot_system(self, plot_style='all', colorby=None, 
        cmap = 'viridis', 
        radius=10, 
        opacity=1.0,
        hide_zero=False,
        color = '#ff7f00'):
        
        if plot_style == 'all':
            visualize.plot_simple(self)
        
        elif plot_style == 'selection':
            visualize.plot_by_selection(self, radius=radius, opacity=opacity)
        
        elif plot_style == 'continuous_property':
            visualize.plot_by_property(sys, colorby, 
                cmap = cmap, 
                radius=radius, 
                opacity=opacity,
                hide_zero=hide_zero)
            
        elif plot_style == 'boolean_property':
            visualize.plot_by_boolean(sys, colorby, 
                color = color, 
                radius=radius, 
                opacity=opacity,
                hide_zero=hide_zero)

    @property
    def material(self):
        if self._material is None:
            self._material = self.graph.value(self.sample, CMSO.hasMaterial)
        return self._material

    @material.setter
    def material(self, value):
        self._material = value

    def duplicate(self, only_essential=False):
        new_system = System()
        if only_essential:
            n_dict = {'positions': copy.deepcopy(self.atoms.positions),
                    'species': copy.deepcopy(self.atoms.species),
                    'types': copy.deepcopy(self.atoms.types),}
        else:
            n_dict = {key: copy.deepcopy(val)[:self.natoms] for key, val in self.atoms.items()}
            new_system.label = self.label
            new_system._name = self._name
        
        atoms = Atoms()
        atoms.from_dict(n_dict)
        atoms._lattice = self.atoms._lattice
        atoms._lattice_constant = self.atoms._lattice_constant
        new_system._structure_dict = copy.deepcopy(self._structure_dict)
        new_system.box = self.box
        new_system.atoms = atoms
        new_system.graph = self.graph
        new_system.sample = None

        return new_system

    def repeat(self, repetitions):
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
        new_system = self.duplicate()
        new_system = operations.repeat(new_system, repetitions)
        if new_system._structure_dict is None:
            new_system._structure_dict = {}
        new_system._structure_dict["repetitions"] = repetitions
        new_system.to_graph()
        new_system.copy_defects(self.sample)
        return new_system
    
    def delete(self, ids=None, indices=None, condition=None, selection=False, copy_structure=False):
        """
        Delete atoms from the structure.

        Parameters
        ----------
        ids : list, optional
            A list of atom IDs to delete. Default is None.
        indices : list, optional
            A list of atom indices to delete. Default is None.
        condition : str, optional
            A condition to select atoms to delete. Default is None.
        selection : bool, optional
            If True, delete atoms based on the current selection. Default is False.
        copy_structure: bool, optional
            If True, a copy of the structure will be returned. Default is False.

        Returns
        -------
        None

        Notes
        -----
        Deletes atoms from the structure based on the provided IDs, indices, condition, or selection.
        If the structure has a graph associated with it, the graph will be updated accordingly.
        """
        if copy_structure:
            sys = self.duplicate()
            #and add this new structure to the graph
            sys.to_graph()
            sys.copy_defects(self.sample)
        else:
            sys = self
        
        masks = sys.atoms._generate_bool_list(
            ids=ids, indices=indices, condition=condition, selection=selection
        )
        
        delete_list = [masks[sys.atoms["head"][x]] for x in range(sys.atoms.ntotal)]
        delete_ids = [x for x in range(sys.atoms.ntotal) if delete_list[x]]
        actual_natoms = sys.natoms
        sys.atoms._delete_atoms(delete_ids)

        if sys.graph is not None:
            # first annotate graph
            val = len([x for x in masks if x])
            c = val / actual_natoms
            sys.add_vacancy(c, number=val)
            # now we need to re-add atoms, so at to remove
            sys.graph.remove((sys.sample, CMSO.hasNumberOfAtoms, None))
            sys.graph.add(
                (
                    sys.sample,
                    CMSO.hasNumberOfAtoms,
                    Literal(actual_natoms - val, datatype=XSD.integer),
                )
            )
            # revamp composition
            # remove existing chem composution

            chemical_species = sys.graph.value(sys.sample, CMSO.hasSpecies)
            # start by cleanly removing elements
            for s in sys.graph.triples((chemical_species, CMSO.hasElement, None)):
                element = s[2]
                sys.graph.remove((element, None, None))
            sys.graph.remove((chemical_species, None, None))
            sys.graph.remove((sys.sample, CMSO.hasSpecies, None))

            # now recalculate and add it again
            composition = sys.schema.material.element_ratio()
            valid = False
            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    valid = True
                    break

            if valid:
                chemical_species = sys.graph.create_node(
                    f"{sys._name}_ChemicalSpecies", CMSO.ChemicalSpecies
                )
                sys.graph.add((sys.sample, CMSO.hasSpecies, chemical_species))

                for e, r in composition.items():
                    if e in element_indetifiers.keys():
                        element = sys.graph.create_node(
                            element_indetifiers[e], CMSO.ChemicalElement
                        )
                        sys.graph.add((chemical_species, CMSO.hasElement, element))
                        sys.graph.add(
                            (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                        )
                        sys.graph.add(
                            (
                                element,
                                CMSO.hasElementRatio,
                                Literal(r, datatype=XSD.float),
                            )
                        )

            # we also have to read in file and clean it up
            filepath = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasPath
            ).toPython()
            position_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasIdentifier
            ).toPython()
            species_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Species"), CMSO.hasIdentifier
            ).toPython()

            # clean up items
            datadict = {
                position_identifier: {
                    "value": sys.schema.atom_attribute.position(),
                    "label": "position",
                },
                species_identifier: {
                    "value": sys.schema.atom_attribute.species(),
                    "label": "species",
                },
            }
            outfile = os.path.join(
                sys.graph.structure_store, str(sys._name).split(":")[-1]
            )
            json_io.write_file(outfile, datadict)

            #write mapping for the operation
            if self.sample.toPython() != sys.sample.toPython():
                activity = self.graph.create_node(f"activity:{uuid.uuid4()}", ASMO.DeleteAtom)
                sys.graph.add((sys.sample, PROV.wasDerivedFrom, self.sample))
                sys.graph.add((sys.sample, PROV.wasGeneratedBy, activity))

        return sys
        
    def add_property_mappings(self, output_property, mapping_quantity=None):
        if self.graph is None:
            return
        if not isinstance(output_property, Property):
            return
        
        #if the property is a directly calculated value
        parent_samples = list([x[0] for x in self.graph.triples((None, CMSO.hasCalculatedProperty, output_property._parent))])
        if len(parent_samples)>0:
            for parent_sample in parent_samples:
                self.graph.add((self.sample, PROV.wasDerivedFrom, parent_sample))
        else:
            #this is quantity that is derived -> for example volume/3 -> it has only sample parent, but no direct connection
            if output_property._sample_parent is not None:
                self.graph.add((self.sample, PROV.wasDerivedFrom, output_property._sample_parent))

        if mapping_quantity=='lattice_constant':
            #add lattice constant mapping
            material = self.graph.value(self.sample, CMSO.hasMaterial)
            crystal_structure = self.graph.value(material, CMSO.hasStructure)
            unit_cell = self.graph.value(crystal_structure, CMSO.hasUnitCell)
            lattice_parameter = self.graph.value(unit_cell, CMSO.hasLatticeParameter)

            #also get activity
            activity = self.graph.value(output_property._parent, ASMO.wasCalculatedBy)
            self.graph.add((lattice_parameter, UNSAFEASMO.wasCalculatedBy, activity))


    def add_vacancy(self, concentration, number=None):
        """
        Add Vacancy details which will be annotated by PODO

        Parameters
        ----------
        concentration: float
            vacancy concentration, value should be between 0-1

        number: int
            Number of atoms that were deleted, optional

        Returns
        -------
        None
        """
        if self.graph is None:
            return

        vacancy = self.graph.create_node(f"{self._name}_Vacancy", PODO.Vacancy)
        self.graph.add((self.material, CMSO.hasDefect, vacancy))
        self.graph.add(
            (
                self.sample,
                PODO.hasVacancyConcentration,
                Literal(concentration, datatype=XSD.float),
            )
        )
        if number is not None:
            self.graph.add(
                (
                    self.sample,
                    PODO.hasNumberOfVacancies,
                    Literal(number, datatype=XSD.integer),
                )
            )

    def substitute_atoms(
        self,
        substitution_element,
        ids=None,
        indices=None,
        condition=None,
        selection=False,
        copy_structure=False,
    ):
        """
        Substitute atoms in the structure with a given element.

        Parameters
        ----------
        substitution_element : str
            The element to substitute the atoms with.
        ids : list, optional
            A list of atom IDs to consider for substitution. Defaults to None.
        indices : list, optional
            A list of atom indices to consider for substitution. Defaults to None.
        condition : callable, optional
            A callable that takes an atom as input and returns a boolean indicating whether the atom should be considered for substitution. Defaults to None.
        selection : bool, optional
            If True, only selected atoms will be considered for substitution. Defaults to False.
        copy_structure: bool, optional
            If True, a copy of the structure will be returned. Defaults to False.

        Returns
        -------
        None

        Notes
        -----
        - This method substitutes atoms in the structure with a given element.
        - The substitution is performed based on the provided IDs, indices, condition, and selection parameters.
        - The substituted atoms will have their species and types updated accordingly.
        - If the graph is not None, the method also operates on the graph by removing existing elements and adding new ones based on the composition of the substituted atoms.
        - The method also cleans up items in the file associated with the graph.

        Examples
        --------
        # Substitute selected atoms with nitrogen
        structure.substitute_atoms("N", ids=[1, 3, 5])
        """
        if copy_structure:
            sys = self.duplicate()
            #and add this new structure to the graph
            sys.to_graph()
            sys.copy_defects(self.sample)
        else:
            sys = self
        
        masks = sys.atoms._generate_bool_list(
            ids=ids, indices=indices, condition=condition, selection=selection
        )
        delete_list = [masks[sys.atoms["head"][x]] for x in range(sys.atoms.ntotal)]
        delete_ids = [x for x in range(sys.atoms.ntotal) if delete_list[x]]
        type_dict = sys.atoms._type_dict
        rtype_dict = {val: key for key, val in type_dict.items()}
        if substitution_element in rtype_dict.keys():
            atomtype = rtype_dict[substitution_element]
            maxtype = atomtype
        else:
            maxtype = max(sys.atoms["types"]) + 1

        for x in delete_ids:
            sys.atoms["species"][x] = substitution_element
            sys.atoms["types"][x] = maxtype
        #impurity metrics
        no_of_impurities = len(delete_ids)
        conc_of_impurities = no_of_impurities/sys.natoms

        # operate on the graph
        if sys.graph is not None:
            chemical_species = sys.graph.value(sys.sample, CMSO.hasSpecies)
            # start by cleanly removing elements
            for s in sys.graph.triples((chemical_species, CMSO.hasElement, None)):
                element = s[2]
                sys.graph.remove((element, None, None))
            sys.graph.remove((chemical_species, None, None))
            sys.graph.remove((sys.sample, CMSO.hasSpecies, None))

            composition = sys.schema.material.element_ratio()
            valid = False
            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    valid = True
                    break

            if valid:
                chemical_species = sys.graph.create_node(
                    f"{sys._name}_ChemicalSpecies", CMSO.ChemicalSpecies
                )
                sys.graph.add((sys.sample, CMSO.hasSpecies, chemical_species))

                for e, r in composition.items():
                    if e in element_indetifiers.keys():
                        element = sys.graph.create_node(
                            element_indetifiers[e], CMSO.ChemicalElement
                        )
                        sys.graph.add((chemical_species, CMSO.hasElement, element))
                        sys.graph.add(
                            (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                        )
                        sys.graph.add(
                            (
                                element,
                                CMSO.hasElementRatio,
                                Literal(r, datatype=XSD.float),
                            )
                        )

            # we also have to read in file and clean it up
            filepath = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasPath
            ).toPython()
            position_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasIdentifier
            ).toPython()
            species_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Species"), CMSO.hasIdentifier
            ).toPython()

            # clean up items
            datadict = {
                position_identifier: {
                    "value": sys.schema.atom_attribute.position(),
                    "label": "position",
                },
                species_identifier: {
                    "value": sys.schema.atom_attribute.species(),
                    "label": "species",
                },
            }
            outfile = os.path.join(
                sys.graph.structure_store, str(sys._name).split(":")[-1]
            )
            json_io.write_file(outfile, datadict)
            sys.add_triples_for_substitutional_impurities(conc_of_impurities, no_of_impurities=no_of_impurities)

            #write mapping for the operation
            if self.sample.toPython() != sys.sample.toPython():
                activity = self.graph.create_node(f"activity:{uuid.uuid4()}", ASMO.SubstituteAtom)
                sys.graph.add((sys.sample, PROV.wasDerivedFrom, self.sample))
                sys.graph.add((sys.sample, PROV.wasGeneratedBy, activity))

        return sys

    def add_triples_for_substitutional_impurities(self, conc_of_impurities, no_of_impurities=None):
        defect = self.graph.create_node(f"{self._name}_SubstitutionalImpurity", PODO.SubstitutionalImpurity)
        self.graph.add((self.material, CMSO.hasDefect, defect))
        self.graph.add((self.sample, PODO.hasImpurityConcentration, Literal(conc_of_impurities, datatype=XSD.float)))
        if no_of_impurities is not None:
            self.graph.add((self.sample, PODO.hasNumberOfImpurityAtoms, Literal(no_of_impurities, datatype=XSD.integer)))

    def add_interstitial_impurities(
        self, element, void_type="tetrahedral",
        lattice_constant=None,
        threshold=0.01,
        copy_structure=False,
    ):
        """
        Add interstitial impurities to the System

        Parameters
        ----------
        element: string or list
            Chemical symbol of the elements/elements to be added
            `element = 'Al'` will add one interstitial while `element = ['Al', 'Al']` or `element = ['Al', 'Li']` will add
            two impurities

        void_type: string
            type of void to be added.  {`tetrahedral`, `octahedral`}

        lattice_constant: float, optional
            lattice constant of the system. Required only for octahedral voids

        threshold: float, optional
            threshold for the distance from the lattice constant for octahedral voids to account for fluctuations in atomic positions
        
        copy_structure: bool, optional
            If True, a copy of the structure will be returned. Defaults to False.

        Returns
        -------
        System:
            system with the added impurities

        Notes
        -----
        The validity of the void positions are not checked! This means that temperature, presence of vacancies or other
        interstitials could affect the addition.
        """
        if None in self.atoms.species:
            raise ValueError("Assign species!")

        sys = self.duplicate()

        if void_type == "tetrahedral":
            element = np.atleast_1d(element)
            self.find.neighbors(method="voronoi", cutoff=0.1)
            verts = self.unique_vertices
            randindex = np.random.randint(0, len(verts), len(element))
            randpos = np.array(verts)[randindex]


        elif void_type == "octahedral":
            if lattice_constant is None:
                if "lattice_constant" in self.lattice_properties.keys():
                    lattice_constant = self.lattice_properties["lattice_constant"]
                else:
                    raise ValueError(
                        "lattice constant is needed for octahedral voids, please provide"
                    )

            cutoff = lattice_constant + threshold * 2
            self.find.neighbors(method="cutoff", cutoff=cutoff)
            octa_pos = []
            for count, dist in enumerate(self.atoms.neighbors.distance):
                diffs = np.abs(np.array(dist) - lattice_constant)
                # print(diffs)
                indices = np.where(diffs < 1e-2)[0]
                # index_neighbor = np.array(self.atoms["neighbors"][count])[indices]
                # real_indices = np.array(self.atoms.neighbors.index[count])[indices]
                # create a dict
                # index_dict = {str(x):y for x,y in zip(real_indices, ghost_indices)}
                vector = np.array(self.atoms["diff"][count])[indices]
                vector = self.atoms.positions[count] + vector / 2
                for vect in vector:
                    vect = self.modify.remap_position_to_box(vect)
                    # print(vect)
                    octa_pos.append(vect)

            octa_pos = np.unique(octa_pos, axis=0)
            randindex = np.random.randint(0, len(octa_pos), len(element))
            randpos = octa_pos[randindex]

            if not len(randpos) == len(element):
                raise ValueError("not enough octahedral positions found!")

        else:
            raise ValueError("void_type can only be tetrahedral/octahedral")

        # create new system with the atoms added
        no_of_impurities = len(randpos)
        conc_of_impurities = no_of_impurities/self.natoms

        if copy_structure:
            #sys = self.duplicate()
            sys = System(source=sys.add_atoms({"positions": randpos, "species": element}))
            sys.graph = self.graph        
            sys.to_graph()
            sys.copy_defects(self.sample)
        else:
            #sys = self.duplicate()
            sys = System(source=self.add_atoms({"positions": randpos, "species": element}))
            sys.graph = self.graph
            sys.sample = self.sample

        # now we have to verify the triples correctly and add them in
        if sys.graph is not None:
            sys.graph.remove((sys.sample, CMSO.hasNumberOfAtoms, None))
            sys.graph.add(
                (
                    sys.sample,
                    CMSO.hasNumberOfAtoms,
                    Literal(sys.natoms, datatype=XSD.integer),
                )
            )
            # revamp composition
            # remove existing chem composution
            chemical_species = sys.graph.value(sys.sample, CMSO.hasSpecies)
            # start by cleanly removing elements
            for s in sys.graph.triples((chemical_species, CMSO.hasElement, None)):
                element = s[2]
                sys.graph.remove((element, None, None))
            sys.graph.remove((chemical_species, None, None))
            sys.graph.remove((sys.sample, CMSO.hasSpecies, None))

            composition = sys.schema.material.element_ratio()
            valid = False
            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    valid = True
                    break

            if valid:
                chemical_species = sys.graph.create_node(
                    f"{sys._name}_ChemicalSpecies", CMSO.ChemicalSpecies
                )
                sys.graph.add((sys.sample, CMSO.hasSpecies, chemical_species))

                for e, r in composition.items():
                    if e in element_indetifiers.keys():
                        element = sys.graph.create_node(
                            element_indetifiers[e], CMSO.ChemicalElement
                        )
                        sys.graph.add((chemical_species, CMSO.hasElement, element))
                        sys.graph.add(
                            (element, CMSO.hasChemicalSymbol, Literal(e, datatype=XSD.string))
                        )
                        sys.graph.add(
                            (
                                element,
                                CMSO.hasElementRatio,
                                Literal(r, datatype=XSD.float),
                            )
                        )

            # we also have to read in file and clean it up
            filepath = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasPath
            ).toPython()
            position_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Position"), CMSO.hasIdentifier
            ).toPython()
            species_identifier = sys.graph.value(
                URIRef(f"{sys.sample}_Species"), CMSO.hasIdentifier
            ).toPython()

            # clean up items
            datadict = {
                position_identifier: {
                    "value": sys.schema.atom_attribute.position(),
                    "label": "position",
                },
                species_identifier: {
                    "value": sys.schema.atom_attribute.species(),
                    "label": "species",
                },
            }
            outfile = os.path.join(
                sys.graph.structure_store, str(sys._name).split(":")[-1]
            )
            json_io.write_file(outfile, datadict)

            sys.add_triples_for_interstitial_impurities(conc_of_impurities, no_of_impurities=no_of_impurities, label=void_type)
            #write mapping for the operation
            if self.sample.toPython() != sys.sample.toPython():
                activity = self.graph.create_node(f"activity:{uuid.uuid4()}", ASMO.AddAtom)
                sys.graph.add((sys.sample, PROV.wasDerivedFrom, self.sample))
                sys.graph.add((sys.sample, PROV.wasGeneratedBy, activity))
        return sys

    def add_triples_for_interstitial_impurities(self, conc_of_impurities, no_of_impurities=None, label=None):
        if label is not None:
            defect = self.graph.create_node(f"{self._name}_InterstitialImpurity", PODO.InterstitialImpurity, label=label)
        else:
            defect = self.graph.create_node(f"{self._name}_InterstitialImpurity", PODO.InterstitialImpurity)
        self.graph.add((self.material, CMSO.hasDefect, defect))
        self.graph.add((self.sample, PODO.hasImpurityConcentration, Literal(conc_of_impurities, datatype=XSD.float)))
        if no_of_impurities is not None:
            self.graph.add((self.sample, PODO.hasNumberOfImpurityAtoms, Literal(no_of_impurities, datatype=XSD.integer)))


    def __delitem__(self, val):
        """
        Delete item(s) from the structure.

        Parameters
        ----------
        val : int or list of int
            The index(es) of the item(s) to be deleted.

        Notes
        -----
        If `val` is an integer, it is converted to a list with a single element.
        The graph is then updated accordingly based on the deleted indices.

        """
        if isinstance(val, int):
            val = [val]
        # now the graph has to be updated accordingly
        self.delete(indices=list(val))

    def write_poscar_id(self, outfile):
        lines = []
        with open(outfile, "r") as fin:
            for line in fin:
                lines.append(line)
        lines[0] = self.sample.toPython() + "\n"
        with open(outfile, "w") as fout:
            for line in lines:
                fout.write(line)

    def to_file(self, filename=None, 
                format='lammps-dump', 
                customkeys=None, customvals=None,
                compressed=False, timestep=0, species=None,  add_sample_id=True,
                copy_from=None, pseudo_files=None):
        """
        Write the structure to a file in the specified format.

        Parameters
        ----------
        outfile : str
            The path to the output file.
        format : str, optional
            The format of the output file. Defaults to 'lammps-dump'.
        customkeys : list, optional
            A list of custom keys to include in the output file. Defaults to None.
            Only valid if format is 'lammps-dump'.
        customvals : list, optional
            A list of custom values corresponding to the custom keys. Defaults to None.
            Only valid if format is 'lammps-dump'.
        compressed : bool, optional
            Whether to compress the output file. Defaults to False.
        timestep : int, optional
            The timestep value to include in the output file. Defaults to 0.
            Only valid if format is 'lammps-dump'.
        species : list, optional
            A list of species to include in the output file. Defaults to None.
            Only valid for ASE, if species is not specified.
        add_sample_id : bool, optional
            Whether to add a sample ID to the output file. Defaults to True.
            Only valid for poscar and quantum-espresso formats.
        copy_from : str, optional
            If provided, input options for quantum-espresso format will be copied from
            the given file. Structure specific information will be replaced.
            Note that the validity of input file is not checked.
        pseudo_files : list, optional
            if provided, add the pseudopotential filenames to file.
            Should be in alphabetical order of chemical species symbols.        

        Returns
        -------
        None
        """
        if filename is None:
            outfile = f"structure.out"
        else:
            outfile = filename

        if format == "ase":
            asesys = convert_snap(self)
            if self.sample is not None:
                asesys.info["sample_id"] = self.sample
            return asesys
        
        elif format == "pyiron":
            from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
            asesys = convert_snap(self)
            pyironsys = ase_to_pyiron(asesys)
            if self.sample is not None:
                pyironsys.info["sample_id"] = self.sample            
            return pyironsys

        elif format == "poscar":
            asesys = convert_snap(self)
            write(outfile, asesys, format="vasp")
            if add_sample_id and (self.sample is not None):
                self.write_poscar_id(outfile)
        
        elif format == "lammps-dump":
            inputmethods.to_file(self, outfile, format='lammps-dump', customkeys=customkeys, customvals=customvals,
                compressed=compressed, timestep=timestep, species=species)
        
        elif format == "lammps-data":
            asesys = convert_snap(self)
            write(outfile, asesys, format='lammps-data', atom_style='atomic')
        
        elif format == "quantum-espresso":
            aio.write_espresso(self, filename, copy_from=copy_from, pseudo_files=pseudo_files)
        
        else:
            asesys = convert_snap(self)
            write(outfile, asesys, format=format)

    def to_graph(self):
        """
        Converts the structure object to a graph representation.

        Returns
        -------
        None
        """
        if self.graph is None:
            return

        self._generate_name()
        self._add_sample()
        self._add_material()
        self._add_chemical_composition()
        self._add_simulation_cell()
        self._add_simulation_cell_properties()
        self._add_crystal_structure()
        self._add_atoms()

    def _generate_name(self, name_index=None):
        if self.names:
            if name_index is None:
                name_index = self.graph.n_samples + 1
            self._name = f"sample:{name_index}"
        else:
            self._name = f"sample:{str(uuid.uuid4())}"

    def _add_sample(self):
        sample = self.graph.create_node(self._name, CMSO.AtomicScaleSample, label=self.label)
        self.sample = sample

    def _add_material(self):
        """
        Add a CMSO Material object

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        material = self.graph.create_node(
            f"{self._name}_Material", CMSO.CrystallineMaterial
        )
        self.graph.add((self.sample, CMSO.hasMaterial, material))
        self.material = material

    def _add_chemical_composition(self):
        """
        Add chemical composition

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        composition = self.schema.material.element_ratio()
        valid = False
        for e, r in composition.items():
            if e in element_indetifiers.keys():
                valid = True
                break

        if valid:
            chemical_species = self.graph.create_node(
                f"{self._name}_ChemicalSpecies", CMSO.ChemicalSpecies
            )
            self.graph.add((self.sample, CMSO.hasSpecies, chemical_species))

            for e, r in composition.items():
                if e in element_indetifiers.keys():
                    element = self.graph.create_node(
                        element_indetifiers[e], CMSO.ChemicalElement
                    )
                    self.graph.add((chemical_species, CMSO.hasElement, element))
                    self.graph.add(
                        (
                            element,
                            CMSO.hasChemicalSymbol,
                            Literal(e, datatype=XSD.string),
                        )
                    )
                    self.graph.add(
                        (element, CMSO.hasElementRatio, Literal(r, datatype=XSD.float))
                    )

    def _add_simulation_cell(self):
        """
        Add a CMSO SimulationCell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        simulation_cell = self.graph.create_node(
            f"{self._name}_SimulationCell", CMSO.SimulationCell
        )
        self.graph.add((self.sample, CMSO.hasSimulationCell, simulation_cell))
        volume = self.graph.create_node(
            f"{self._name}_Volume", UNSAFEASMO.Volume, label="SimulationCellVolume"
        )
        self.graph.add((simulation_cell, UNSAFECMSO.hasVolume, volume))
        self.graph.add(
            (
                volume,
                UNSAFEASMO.hasValue,
                Literal(
                    np.round(self.schema.simulation_cell.volume(), decimals=2),
                    datatype=XSD.float,
                ),
            )
        )
        self.graph.add(
            (
                volume,
                UNSAFEASMO.hasUnit,
                URIRef(f"http://qudt.org/vocab/unit/ANGSTROM3"),
            )
        )
        self.graph.add(
            (
                self.sample,
                CMSO.hasNumberOfAtoms,
                Literal(
                    self.schema.simulation_cell.number_of_atoms(), datatype=XSD.integer
                ),
            )
        )
        
        repetitions = self.schema.simulation_cell.repetitions()
        self.graph.add((simulation_cell, CMSO.hasRepetition_x, Literal(repetitions[0], datatype=XSD.integer)))
        self.graph.add((simulation_cell, CMSO.hasRepetition_y, Literal(repetitions[1], datatype=XSD.integer)))
        self.graph.add((simulation_cell, CMSO.hasRepetition_z, Literal(repetitions[2], datatype=XSD.integer)))
        self.simulation_cell = simulation_cell

    def _add_simulation_cell_properties(self):
        """
        Add a CMSO SimulationCell properties such as SimulationCellLength,
        and Vectors.

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        simulation_cell_length = self.graph.create_node(
            f"{self._name}_SimulationCellLength", CMSO.SimulationCellLength
        )
        self.graph.add((self.simulation_cell, CMSO.hasLength, simulation_cell_length))
        data = self.schema.simulation_cell.length()
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_x,
                Literal(data[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_y,
                Literal(data[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_length,
                CMSO.hasLength_z,
                Literal(data[2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_01 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_1", CMSO.SimulationCellVector
        )
        data = self.schema.simulation_cell.vector()
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_01)
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_x,
                Literal(data[0][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_y,
                Literal(data[0][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_01,
                CMSO.hasComponent_z,
                Literal(data[0][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_02 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_2", CMSO.SimulationCellVector
        )
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_02)
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_x,
                Literal(data[1][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_y,
                Literal(data[1][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_02,
                CMSO.hasComponent_z,
                Literal(data[1][2], datatype=XSD.float),
            )
        )

        simulation_cell_vector_03 = self.graph.create_node(
            f"{self._name}_SimulationCellVector_3", CMSO.SimulationCellVector
        )
        self.graph.add(
            (self.simulation_cell, CMSO.hasVector, simulation_cell_vector_03)
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_x,
                Literal(data[2][0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_y,
                Literal(data[2][1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_vector_03,
                CMSO.hasComponent_z,
                Literal(data[2][2], datatype=XSD.float),
            )
        )

        simulation_cell_angle = self.graph.create_node(
            f"{self._name}_SimulationCellAngle", CMSO.SimulationCellAngle
        )
        data = self.schema.simulation_cell.angle()
        self.graph.add((self.simulation_cell, CMSO.hasAngle, simulation_cell_angle))
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_alpha,
                Literal(data[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_beta,
                Literal(data[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                simulation_cell_angle,
                CMSO.hasAngle_gamma,
                Literal(data[2], datatype=XSD.float),
            )
        )

    def _add_crystal_structure(self, targets=None):
        """
        Add a CMSO Crystal Structure

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        if targets is None:
            targets = [
                self.schema.material.crystal_structure.name(),
                self.schema.material.crystal_structure.spacegroup_symbol(),
                self.schema.material.crystal_structure.spacegroup_number(),
                self.schema.material.crystal_structure.unit_cell.bravais_lattice(),
                self.schema.material.crystal_structure.unit_cell.lattice_parameter(),
                self.schema.material.crystal_structure.unit_cell.angle(),
            ]

        #fix for lattice angle of HCP
        if targets[0] == 'hcp':
            targets[5] = [90.0, 90.0, 120.0]

        valid = self.graph._is_valid(targets)

        if valid:
            crystal_structure = self.graph.create_node(
                f"{self._name}_CrystalStructure", CMSO.CrystalStructure
            )
            self.graph.add((self.material, CMSO.hasStructure, crystal_structure))
            self.graph.add(
                (
                    crystal_structure,
                    CMSO.hasAltName,
                    Literal(targets[0], datatype=XSD.string),
                )
            )
            self.crystal_structure = crystal_structure

            if targets[1] is not None:
                self._add_space_group(targets[1], targets[2])

            # now see if unit cell needs to be added
            valid = self.graph._is_valid(targets[3:])
            if valid:
                self._add_unit_cell()
                if targets[3] is not None:
                    self._add_bravais_lattice(targets[3])
                if targets[4] is not None:
                    self._add_lattice_properties(targets[4], targets[5])

    def _add_space_group(self, spacegroup_symbol, spacegroup_number):
        """
        Add a CMSO Space Group

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        #space_group = URIRef(f"{self._name}_SpaceGroup")
        #self.graph.add((self.crystal_structure, CMSO.hasSpaceGroup, space_group))
        self.graph.add(
            (
                self.crystal_structure,
                CMSO.hasSpaceGroupSymbol,
                Literal(spacegroup_symbol, datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                self.crystal_structure,
                CMSO.hasSpaceGroupNumber,
                Literal(spacegroup_number, datatype=XSD.integer),
            )
        )

    def _add_unit_cell(self):
        """
        Add a CMSO Unit Cell

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """

        unit_cell = self.graph.create_node(f"{self._name}_UnitCell", CMSO.UnitCell)
        self.graph.add((self.crystal_structure, CMSO.hasUnitCell, unit_cell))
        self.unit_cell = unit_cell

    def _add_bravais_lattice(self, bv):
        """
        Add a Bravais lattice to the unit cell.

        Parameters:
            bv (str): The URI of the Bravais lattice.

        Returns:
            None
        """
        bv = URIRef(bv)
        self.graph.add(
            (
                self.unit_cell,
                CMSO.hasBravaisLattice,
                bv,
            )
        )

    def _add_lattice_properties(self, lattice_parameter_value, lattice_angle_value):
        """
        Add CMSO lattice properties such as Lattice Parameter,
        and its lengths and angles.

        Parameters
        ----------
        name
            if provided, the name will be used instead of random identifier

        Returns
        -------
        """
        lattice_parameter = self.graph.create_node(
            f"{self._name}_LatticeParameter", CMSO.LatticeParameter
        )
        self.graph.add((self.unit_cell, CMSO.hasLatticeParameter, lattice_parameter))
        self.graph.add(
            (
                lattice_parameter,
                CMSO.hasLength_x,
                Literal(lattice_parameter_value[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_parameter,
                CMSO.hasLength_y,
                Literal(lattice_parameter_value[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_parameter,
                CMSO.hasLength_z,
                Literal(lattice_parameter_value[2], datatype=XSD.float),
            )
        )

        lattice_angle = self.graph.create_node(
            f"{self._name}_LatticeAngle", CMSO.LatticeAngle
        )
        self.graph.add((self.unit_cell, CMSO.hasAngle, lattice_angle))
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_alpha,
                Literal(lattice_angle_value[0], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_beta,
                Literal(lattice_angle_value[1], datatype=XSD.float),
            )
        )
        self.graph.add(
            (
                lattice_angle,
                CMSO.hasAngle_gamma,
                Literal(lattice_angle_value[2], datatype=XSD.float),
            )
        )

    def _save_atom_attributes(self, position_identifier, species_identifier):
        """
        Save the atom attributes to a file.

        Parameters
        ----------
        position_identifier : str
            The identifier for the position attribute.
        species_identifier : str
            The identifier for the species attribute.

        Returns
        -------
        str
            The relative path to the saved file.

        Notes
        -----
        This method saves the atom attributes to a file in the file-based store system.
        The attributes are stored in a dictionary with the position identifier and species identifier as keys.
        The dictionary is then written to a JSON file using the `json_io.write_file` function.
        The file is saved in the structure store directory with the name of the structure as the filename.
        The method returns the relative path to the saved file.
        """
        datadict = {
            position_identifier: {
                "value": self.schema.atom_attribute.position(),
                "label": "position",
            },
            species_identifier: {
                "value": self.schema.atom_attribute.species(),
                "label": "species",
            },
        }
        outfile = os.path.join(
            self.graph.structure_store, str(self._name).split(":")[-1]
        )
        json_io.write_file(outfile, datadict)
        return os.path.relpath(outfile + ".json")

    def _add_atoms(self):
        """
        Add Atoms including their species and positions

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Note that for the moment, we will dump the structures in a given folder,
        maybe this could be input from the Job class directly
        """
        # now we write out file
        position_identifier = str(uuid.uuid4())
        species_identifier = str(uuid.uuid4())

        outfile = self._save_atom_attributes(position_identifier, species_identifier)

        if "positions" in self.atoms.keys():
            position = self.graph.create_node(
                f"{self._name}_Position", CMSO.AtomAttribute
            )
            self.graph.add(
                (
                    self.sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    position,
                )
            )
            self.graph.add(
                (position, CMSO.hasName, Literal("Position", datatype=XSD.string))
            )
            self.graph.add(
                (
                    position,
                    CMSO.hasIdentifier,
                    Literal(position_identifier, datatype=XSD.string),
                )
            )
            self.graph.add(
                (position, CMSO.hasPath, Literal(outfile, datatype=XSD.string))
            )

        if "species" in self.atoms.keys():
            species = self.graph.create_node(
                f"{self._name}_Species", CMSO.AtomAttribute
            )
            self.graph.add(
                (
                    self.sample,
                    Namespace("http://purls.helmholtz-metadaten.de/cmso/").hasAttribute,
                    species,
                )
            )
            self.graph.add(
                (species, CMSO.hasName, Literal("Species", datatype=XSD.string))
            )
            self.graph.add(
                (
                    species,
                    CMSO.hasIdentifier,
                    Literal(species_identifier, datatype=XSD.string),
                )
            )
            self.graph.add(
                (species, CMSO.hasPath, Literal(outfile, datatype=XSD.string))
            )

        # if "velocities" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Velocity'
        #    velocity = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, velocity))
        #    self.add((velocity, RDF.type, CMSO.AtomAttribute))
        #    self.add((velocity, CMSO.hasName, Literal('Velocity', data_type=XSD.string)))
        #    velocity_identifier = uuid.uuid4()
        #    self.add((velocity, CMSO.hasIdentifier, Literal(velocity_identifier, datatype=XSD.string)))

        # if "forces" in self.sys.atoms.keys():
        #    uname = None
        #    if name is not None:
        #        uname = f'{name}_Force'
        #    force = BNode(uname)
        #    self.add((self.sample, CMSO.hasAttribute, force))
        #    self.add((force, RDF.type, CMSO.AtomAttribute))
        #    self.add((force, CMSO.hasName, Literal('Force', data_type=XSD.string)))
        #    force_identifier = uuid.uuid4()
        #    self.add((force, CMSO.hasIdentifier, Literal(force_identifier, datatype=XSD.string)))

    def add_dislocation(self, disl_dict):
        if self.graph is None:
            return
        
        #find what kind of disl is present
        angle_deg = disl_dict['DislocationCharacter']
        if (np.abs(angle_deg-0) < 1E-3) or (np.abs(angle_deg-180) < 1E-3) or (np.abs(angle_deg-360) < 1E-3):
            disl_type = LDO.ScrewDislocation
            disl_name = "ScrewDislocation"
        elif (np.abs(angle_deg-90) < 1E-3) or (np.abs(angle_deg-270) < 1E-3):
            disl_type = LDO.EdgeDislocation
            disl_name = "EdgeDislocation"
        else:
            disl_type = LDO.MixedDislocation
            disl_name = "MixedDislocation"

        line_defect = self.graph.create_node(f"{self._name}_Dislocation", disl_type)
        self.graph.add((self.material, CMSO.hasDefect, line_defect))

        line_direction = self.graph.create_node(f"{self._name}_DislocationLineDirection", LDO.LineDirection)
        self.graph.add((line_direction, CMSO.hasComponent_x, Literal(disl_dict['DislocationLine'][0], datatype=XSD.float)))
        self.graph.add((line_direction, CMSO.hasComponent_y, Literal(disl_dict['DislocationLine'][1], datatype=XSD.float)))
        self.graph.add((line_direction, CMSO.hasComponent_z, Literal(disl_dict['DislocationLine'][2], datatype=XSD.float)))
        self.graph.add((line_defect, LDO.hasLineDirection, line_direction))                

        burgers_vector = self.graph.create_node(f"{self._name}_DislocationBurgersVector", LDO.BurgersVector)
        self.graph.add((burgers_vector, CMSO.hasComponent_x, Literal(disl_dict['BurgersVector'][0], datatype=XSD.float)))
        self.graph.add((burgers_vector, CMSO.hasComponent_y, Literal(disl_dict['BurgersVector'][1], datatype=XSD.float)))
        self.graph.add((burgers_vector, CMSO.hasComponent_z, Literal(disl_dict['BurgersVector'][2], datatype=XSD.float)))
        self.graph.add((line_defect, LDO.hasBurgersVector, burgers_vector))

        if disl_name == "MixedDislocation":
            self.graph.add((line_defect, LDO.hasCharacterAngle, Literal(angle_deg, datatype=XSD.float)))

        slip_direction = self.graph.create_node(f"{self._name}_DislocationSlipDirection", LDO.SlipDirection)
        self.graph.add((slip_direction, CMSO.hasComponent_x, Literal(disl_dict['SlipDirection'][0], datatype=XSD.float)))
        self.graph.add((slip_direction, CMSO.hasComponent_y, Literal(disl_dict['SlipDirection'][1], datatype=XSD.float)))
        self.graph.add((slip_direction, CMSO.hasComponent_z, Literal(disl_dict['SlipDirection'][2], datatype=XSD.float)))
        
        slip_plane = self.graph.create_node(f"{self._name}_DislocationSlipPlane", LDO.SlipPlane)
        normal_vector = self.graph.create_node(f"{self._name}_DislocationNormalVector", LDO.NormalVector)
        self.graph.add((normal_vector, CMSO.hasComponent_x, Literal(disl_dict['SlipPlane'][0], datatype=XSD.float)))
        self.graph.add((normal_vector, CMSO.hasComponent_y, Literal(disl_dict['SlipPlane'][1], datatype=XSD.float)))
        self.graph.add((normal_vector, CMSO.hasComponent_z, Literal(disl_dict['SlipPlane'][2], datatype=XSD.float)))
        self.graph.add((slip_plane, LDO.hasNormalVector, normal_vector))

        slip_system = self.graph.create_node(f"{self._name}_DislocationSlipSystem", LDO.SlipSystem)
        self.graph.add((slip_direction, LDO.belongsToSystem, slip_system))
        self.graph.add((slip_plane, LDO.belongsToSystem, slip_system))
        self.graph.add((line_defect, LDO.movesOn, slip_system))


    def add_stacking_fault(self, sf_dict):
        if self.graph is None:
            return
        plane = " ".join(np.array(sf_dict["plane"]).astype(str))
        displ = " ".join(np.array(sf_dict["displacement"]).astype(str))
        sf = self.graph.create_node(f"{self._name}_StackingFault", PLDO.StackingFault)
        self.graph.add((self.material, CMSO.hasDefect, sf))
        self.graph.add((sf, PLDO.hasSFplane, Literal(plane, datatype=XSD.string)))
        self.graph.add((sf, PLDO.hasDisplacementVector, Literal(displ, datatype=XSD.string)))

    def add_gb(self, gb_dict):
        """
        Add GB details which will be annotated using PLDO

        Parameters
        ----------
        gb_dict : dict
            A dictionary containing details about the grain boundary.
            It should have the following keys:
            - "GBType" (str): The type of grain boundary. Possible values are "Twist", "Tilt", "Symmetric Tilt", and "Mixed".
            - "sigma" (int): The sigma value of the grain boundary.
            - "GBPlane" (str): The plane of the grain boundary.
            - "RotationAxis" (list): The rotation axis of the grain boundary.
            - "MisorientationAngle" (float): The misorientation angle of the grain boundary.

        Returns
        -------
        None

        Notes
        -----
        This method adds grain boundary details to the structure and annotates it using PLDO ontology.
        The grain boundary type, sigma value, GB plane, rotation axis, and misorientation angle are stored as attributes of the grain boundary node in the graph.
        """
        # mark that the structure has a defect
        if self.graph is None:
            return

        if gb_dict["GBType"] is None:
            plane_defect = self.graph.create_node(f"{self._name}_GrainBoundary", PLDO.GrainBoundary)

        elif gb_dict["GBType"] == "Twist":
            plane_defect = self.graph.create_node(
                f"{self._name}_TwistGrainBoundary", PLDO.TwistGrainBoundary
            )

        elif gb_dict["GBType"] == "Tilt":
            plane_defect = self.graph.create_node(
                f"{self._name}_TiltGrainBoundary", PLDO.TiltGrainBoundary
            )

        elif gb_dict["GBType"] == "Symmetric Tilt":
            plane_defect = self.graph.create_node(
                f"{self._name}_SymmetricalTiltGrainBoundary",
                PLDO.SymmetricalTiltGrainBoundary,
            )

        elif gb_dict["GBType"] == "Mixed":
            plane_defect = self.graph.create_node(
                f"{self._name}_MixedGrainBoundary", PLDO.MixedGrainBoundary
            )

        self.graph.add((self.material, CMSO.hasDefect, plane_defect))
        self.graph.add(
            (
                plane_defect,
                PLDO.hasSigmaValue,
                Literal(gb_dict["sigma"], datatype=XSD.integer),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasGBplane,
                Literal(gb_dict["GBPlane"], datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasRotationAxis,
                Literal(gb_dict["RotationAxis"], datatype=XSD.string),
            )
        )
        self.graph.add(
            (
                plane_defect,
                PLDO.hasMisorientationAngle,
                Literal(gb_dict["MisorientationAngle"], datatype=XSD.float),
            )
        )

    
    def rotate(self, rotation_vectors, graph=None, label=None):
        try:
            from atomman.defect.Dislocation import Dislocation
            import atomman as am
            import atomman.unitconvert as uc
        except ImportError:
            raise ImportError("This function requires the atomman package to be installed")
    
        box = am.Box(
            avect=self.box[0],
            bvect=self.box[1],
            cvect=self.box[2],
        )
        
        atoms = am.Atoms(
            atype=self.atoms.types, pos=self.atoms.positions
        )
        
        element = [val for key, val in self.atoms._type_dict.items()]

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
            output_structure.graph = self.graph
        output_structure.atoms._lattice = self.atoms._lattice
        output_structure.atoms._lattice_constant = self.atoms._lattice_constant
        output_structure._structure_dict = self._structure_dict
        if label is not None:
            output_structure.label = label
        else:
            output_structure.label = self.label
        output_structure.to_graph()
        output_structure.copy_defects(self.sample)
        if output_structure.graph is not None:
            self.add_rotation_triples(rotation_vectors, output_structure.sample)
        return output_structure

    def add_rotation_triples(self, rotation_vectors, child_sample_id):
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Rotation)
        self.graph.add((child_sample_id, PROV.wasGeneratedBy, activity))
        self.graph.add((child_sample_id, PROV.wasDerivedFrom, self.sample))

        rot_vector_01 = self.graph.create_node(f"{activity_id}_RotationVector_1", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_01))
        self.graph.add((rot_vector_01, CMSO.hasComponent_x, Literal(rotation_vectors[0][0], datatype=XSD.float),))
        self.graph.add((rot_vector_01, CMSO.hasComponent_y, Literal(rotation_vectors[0][1], datatype=XSD.float),))
        self.graph.add((rot_vector_01, CMSO.hasComponent_z, Literal(rotation_vectors[0][2], datatype=XSD.float),))

        rot_vector_02 = self.graph.create_node(f"{activity_id}_RotationVector_2", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_02))
        self.graph.add((rot_vector_02, CMSO.hasComponent_x, Literal(rotation_vectors[1][0], datatype=XSD.float),))
        self.graph.add((rot_vector_02, CMSO.hasComponent_y, Literal(rotation_vectors[1][1], datatype=XSD.float),))
        self.graph.add((rot_vector_02, CMSO.hasComponent_z, Literal(rotation_vectors[1][2], datatype=XSD.float),))

        rot_vector_03 = self.graph.create_node(f"{activity_id}_RotationVector_3", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, rot_vector_03))
        self.graph.add((rot_vector_03, CMSO.hasComponent_x, Literal(rotation_vectors[2][0], datatype=XSD.float),))
        self.graph.add((rot_vector_03, CMSO.hasComponent_y, Literal(rotation_vectors[2][1], datatype=XSD.float),))
        self.graph.add((rot_vector_03, CMSO.hasComponent_z, Literal(rotation_vectors[2][2], datatype=XSD.float),))

    def _select_by_plane(self, plane, distance, reverse_orientation=False):
        plane_norm = np.linalg.norm(plane)
        selection = []
        for pos in self.atoms.positions:
            dist = np.dot(plane, pos)/plane_norm
            
            if dist < distance:
                selection.append(True)
            else:
                selection.append(False)
        if reverse_orientation:
            selection = np.invert(selection)
        return selection             

    def select_by_plane(self, plane, distance, reverse_orientation=False):
        selection = self._select_by_plane(plane, distance, 
                        reverse_orientation=reverse_orientation)
        self.apply_selection(condition=selection)
        
    def translate(self, translation_vector, 
                        plane=None, distance=None, 
                        reverse_orientation=False, 
                        copy_structure=True,
                        add_triples=True):
        
        if copy_structure:
            sys = self.duplicate()
            #and add this new structure to the graph
            sys.to_graph()
            sys.copy_defects(self.sample)
        else:
            sys = self

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

        if (sys.graph is not None) and add_triples:
            sys.add_translation_triples(translation_vector, plane, distance)
            if self.sample.toPython() != sys.sample.toPython():
                sys.graph.add((sys.sample, PROV.wasDerivedFrom, self.sample))
        return sys
    
    def add_translation_triples(self, translation_vector, plane, distance, ):
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Translation)
        self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

        #now add specifics
        #shear is a vector
        t_vector = self.graph.create_node(f"{activity_id}_TranslationVector", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, t_vector))
        self.graph.add((t_vector, CMSO.hasComponent_x, Literal(translation_vector[0], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_y, Literal(translation_vector[1], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_z, Literal(translation_vector[2], datatype=XSD.float),))

    def shear(self, shear_vector, 
                    plane, 
                    distance, 
                    reverse_orientation=False, 
                    copy_structure=True):
        
        if copy_structure:
            sys = self.duplicate()
            #and add this new structure to the graph
            sys.to_graph()
            sys.copy_defects(self.sample)
        else:
            sys = self

        if not np.dot(shear_vector, plane) == 0:
            raise ValueError("shear vector must be perpendicular to the plane")
        
        sys = sys.translate(shear_vector, plane=plane, distance=distance, 
                            reverse_orientation=reverse_orientation, copy_structure=False,
                            add_triples=False)

        if sys.graph is not None:
            sys.add_shear_triples(shear_vector, plane, distance)
            if self.sample.toPython() != sys.sample.toPython():
                sys.graph.add((sys.sample, PROV.wasDerivedFrom, self.sample))
        return sys
    
    def add_shear_triples(self, translation_vector, plane, distance, ):
        activity_id = f"operation:{uuid.uuid4()}"
        activity = self.graph.create_node(activity_id, ASMO.Shear)
        self.graph.add((self.sample, PROV.wasGeneratedBy, activity))

        #now add specifics
        #shear is a vector
        t_vector = self.graph.create_node(f"{activity_id}_ShearVector", CMSO.Vector)
        self.graph.add((activity, CMSO.hasVector, t_vector))
        self.graph.add((t_vector, CMSO.hasComponent_x, Literal(translation_vector[0], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_y, Literal(translation_vector[1], datatype=XSD.float),))
        self.graph.add((t_vector, CMSO.hasComponent_z, Literal(translation_vector[2], datatype=XSD.float),))

        #if plane is provided, add that as well
        if plane is not None:
            plane = self.graph.create_node(f"{activity_id}_Plane", CMSO.Plane)
            plane_vector = self.graph.create_node(f"{activity_id}_PlaneVector", CMSO.NormalVector)
            self.graph.add((activity, UNSAFECMSO.hasPlane, plane))
            self.graph.add((plane, CMSO.hasNormalVector, plane_vector))
            self.graph.add((plane_vector, CMSO.hasComponent_x, Literal(plane[0], datatype=XSD.float),))
            self.graph.add((plane_vector, CMSO.hasComponent_y, Literal(plane[1], datatype=XSD.float),))
            self.graph.add((plane_vector, CMSO.hasComponent_z, Literal(plane[2], datatype=XSD.float),))
            self.graph.add((plane, CMSO.hasDistanceFromOrigin, Literal(distance, datatype=XSD.float)))

    def copy_defects(self, parent_sample):
        if self.sample is None:
            return
        if parent_sample is None:
            return
        self.graph.copy_defects(self.sample, parent_sample)

    def plot3d(self, *args, **kwargs):
        try:
            from pyiron_atomistics.atomistics.structure.atoms import (
                ase_to_pyiron,
                pyiron_to_ase,
            )
        except ImportError:
            raise ImportError("Please install pyiron_atomistics")
        ase_structure = self.write.ase()
        pyiron_structure = ase_to_pyiron(ase_structure)
        return pyiron_structure.plot3d(*args, **kwargs)
