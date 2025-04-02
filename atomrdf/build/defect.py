import pyscal3.structure_creator as pcs
from atomrdf.system import System
from atomrdf.build.buildutils import _declass
from atomrdf.build.bulk import bulk
import numpy as np
import pyscal3.core as pc
from pyscal3.core import structure_dict, element_dict


def _delete_atom(system, 
    ids=None, 
    indices=None, 
    condition=None, 
    selection=False, 
    copy_structure=False):
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
        sys = system.duplicate()
        #and add this new structure to the graph
        sys.to_graph()
        sys.copy_defects(system.sample)
    else:
        sys = system
    
    masks = sys.atoms._generate_bool_list(
        ids=ids, indices=indices, condition=condition, selection=selection
    )
    
    delete_list = [masks[sys.atoms["head"][x]] for x in range(sys.atoms.ntotal)]
    delete_ids = [x for x in range(sys.atoms.ntotal) if delete_list[x]]
    actual_natoms = sys.natoms
    sys.atoms._delete_atoms(delete_ids)
    vacancy_no = len([x for x in masks if x])
    concentration = vacancy_no / actual_natoms
    sys.add_vacancy(concentration, number=vacancy_no)
    sys.update_system_for_vacancy_creation(vacancy_no,
        actual_natoms,copy_structure=copy_structure)
    return sys

def vacancy(
    ids=None,
    indices=None,
    no_of_vacancies=None,
    structure_object = None,
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
    copy_structure=False,
):
    if all(v is None for v in [ids, indices, no_of_vacancies]):
        raise ValueError("Please provide ids, indices or no_of_vacancies")

    if structure_object is not None:
        #structure is given, so we use that for further work
        sys = structure_object
    else:
        #we try creating a structure, in that case, we call bulk        
        sys = bulk(
            element=element,
            structure=structure,
            lattice_constant=lattice_constant,
            repetitions=repetitions,
            ca_ratio=ca_ratio,
            noise=noise,
            primitive=primitive,
            graph=graph,
            names=names,
            label=label,
        )
    
    if no_of_vacancies is not None:
        indices = np.random.randit(0, sys.natoms, no_of_vacancies)
    
    #get actual natoms
    actual_natoms = sys.natoms

    #we are good to go now
    sys = _delete_atom(sys, 
                ids=ids, 
                indices=indices, 
                condition=None, 
                selection=False, 
                copy_structure=copy_structure)
    return sys


def substitutional(
    system,
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
        sys = system.duplicate()
        #and add this new structure to the graph
        sys.to_graph()
        sys.copy_defects(system.sample)
    else:
        sys = system
    
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
    sys.update_system_for_substitutional_impurities(no_of_impurities,
                    actual_natoms=sys.natoms, copy_structure=copy_structure)
    sys.add_substitutional_impurities(conc_of_impurities,
                                      no_of_impurities=no_of_impurities)

    return sys


def interstitial(
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
    sys.update_system_for_interstitial_impurity(copy_structure=copy_structure)
    sys.add_interstitial_impurities(conc_of_impurities,
                                    no_of_impurities=no_of_impurities)
    return sys

def stacking_fault(
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


def dislocation(
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

def grain_boundary(
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