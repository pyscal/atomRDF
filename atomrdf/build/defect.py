from platform import system
from atomrdf.build.bulk import bulk, _generate_atomic_sample_data
from atomrdf.datamodels.basemodels import DataProperty
import numpy as np
from atomrdf.datamodels.defects.dislocation import Dislocation
from atomrdf.datamodels.defects.stackingfault import StackingFault
from atomrdf.datamodels.defects.pointdefects import (
    Vacancy,
    Interstitial,
    Substitutional,
)
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.build.buildutils import _declass
from pyscal3.grain_boundary import GrainBoundary
from ase import Atoms
import atomrdf.datamodels.workflow.operations as ops


def stacking_fault(
    element,
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
    crystalstructure=None,
    a=None,
    b=None,
    c=None,
    alpha=None,
    covera=None,
    repeat=1,
    graph=None,
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

    a = _declass(a)
    b = _declass(b)
    c = _declass(c)
    alpha = _declass(alpha)
    covera = _declass(covera)

    input_structure, sdict = bulk(
        element,
        crystalstructure=crystalstructure,
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        covera=covera,
        repeat=repeat,
        get_metadata=True,
    )

    ucell = am.load("ase_Atoms", input_structure)
    sf = am.defect.StackingFault(slip_plane, ucell)
    if slip_direction_a is not None:
        sf.a1vect_uvw = slip_direction_a
    if slip_direction_b is not None:
        sf.a2vect_uvw = slip_direction_b
    surfacesystem = sf.surface(
        shift=sf.shifts[0], minwidth=minwidth, even=even, vacuumwidth=vacuum
    )
    if relative_fault_position != 0.5:
        sf.faultpos_rel = relative_fault_position
    faultsystem = sf.fault(a1=displacement_a, a2=displacement_b)

    # get displacements
    displ = am.displacement(surfacesystem, faultsystem)
    aseatoms = faultsystem.dump("ase_Atoms", return_prop=False)
    aseatoms.set_positions(aseatoms.get_positions() + displ)

    if graph is not None:
        data = _generate_atomic_sample_data(aseatoms, sdict, repeat)
        sample = AtomicScaleSample(**data)

        datadict = StackingFault.template()
        datadict["plane"] = slip_plane
        # Convert numpy array to list for pydantic validation
        datadict["displacement"] = displ.flatten().tolist()
        setattr(sample, "stacking_fault", StackingFault(**datadict))
        sample.to_graph(graph)
        aseatoms.info["id"] = sample.id

    return aseatoms


def dislocation(
    element,
    slip_system,
    dislocation_line,
    elastic_constant_dict,
    burgers_vector=None,
    dislocation_type="monopole",
    crystalstructure=None,
    a=None,
    b=None,
    c=None,
    alpha=None,
    covera=None,
    repeat=1,
    graph=None,
    label=None,
    return_atomman_dislocation=False,
):
    """
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
    a = _declass(a)
    b = _declass(b)
    c = _declass(c)
    alpha = _declass(alpha)
    covera = _declass(covera)

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
        raise ValueError("burgers vector should be None, scalar, or of length 3")

    input_structure, sdict = bulk(
        element,
        crystalstructure=crystalstructure,
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        covera=covera,
        repeat=repeat,
        get_metadata=True,
    )

    for key, val in elastic_constant_dict.items():
        elastic_constant_dict[key] = uc.set_in_units(val, "GPa")
    C = am.ElasticConstants(**elastic_constant_dict)

    box = am.Box(
        avect=input_structure.cell[0],
        bvect=input_structure.cell[1],
        cvect=input_structure.cell[2],
    )

    types = [1 for x in range(len(input_structure))]

    atoms = am.Atoms(atype=types, pos=input_structure.get_positions())
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

    aseatoms = disl_system.dump("ase_Atoms", return_prop=False)

    if graph is not None:
        data = _generate_atomic_sample_data(aseatoms, sdict, repeat)
        sample = AtomicScaleSample(**data)

        # now we need to add the dislocation info

        # find dislocation character
        angle = np.dot(dislocation_line, burgers_vector) / (
            np.linalg.norm(dislocation_line) * np.linalg.norm(burgers_vector)
        )
        angle_rad = np.arccos(angle)
        angle_deg = np.degrees(angle_rad)

        if (
            (np.abs(angle_deg - 0) < 1e-3)
            or (np.abs(angle_deg - 180) < 1e-3)
            or (np.abs(angle_deg - 360) < 1e-3)
        ):
            from atomrdf.datamodels.defects.dislocation import (
                ScrewDislocation as DislocationObject,
            )

            disl_name = "screw_dislocation"
        elif (np.abs(angle_deg - 90) < 1e-3) or (np.abs(angle_deg - 270) < 1e-3):
            from atomrdf.datamodels.defects.dislocation import (
                EdgeDislocation as DislocationObject,
            )

            disl_name = "edge_dislocation"
        else:
            from atomrdf.datamodels.defects.dislocation import (
                MixedDislocation as DislocationObject,
            )

            disl_name = "mixed_dislocation"

        disl_dict = DislocationObject.template()
        disl_dict["line_direction"] = dislocation_line
        disl_dict["burgers_vector"] = burgers_vector
        disl_dict["slip_system"]["slip_direction"] = slip_direction
        disl_dict["slip_system"]["slip_plane"]["normal"] = slip_plane
        if disl_name == "mixed_dislocation":
            disl_dict["character_angle"] = angle_deg

        setattr(sample, disl_name, DislocationObject(**disl_dict))
        sample.to_graph(graph)
        aseatoms.info["id"] = sample.id

    if return_atomman_dislocation:
        return aseatoms, disc
    return aseatoms


def grain_boundary(
    element,
    axis,
    sigma,
    gb_plane,
    crystalstructure=None,
    a=None,
    b=None,
    c=None,
    alpha=None,
    covera=None,
    overlap=0.0,
    gap=0.0,
    vacuum=0.0,
    delete_layer="0b0t0b0t",
    tolerance=0.25,
    uc_a=1,
    uc_b=1,
    repeat=None,
    graph=None,
    primitive=False,
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
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
        from aimsgb import GrainBoundary as AIMSGrainBoundary
        from aimsgb import Grain as AIMSGrain
    except ImportError:
        raise ImportError(
            "This function requires the aimsgb and pymatgen packages to be installed"
        )

    a = _declass(a)
    b = _declass(b)
    c = _declass(c)
    alpha = _declass(alpha)
    covera = _declass(covera)

    input_structure, sdict = bulk(
        element,
        crystalstructure=crystalstructure,
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        covera=covera,
        repeat=repeat,
        get_metadata=True,
    )

    pmsys = AseAtomsAdaptor().get_structure(atoms=input_structure)
    grain = AIMSGrain(pmsys.lattice, pmsys.species, pmsys.frac_coords)
    gb = AIMSGrainBoundary(
        axis=axis,
        sigma=sigma,
        plane=gb_plane,
        initial_struct=grain,
        uc_a=uc_a,
        uc_b=uc_b,
    )
    gb_struct = AIMSGrain.stack_grains(
        grain_a=gb.grain_a,
        grain_b=gb.grain_b,
        vacuum=vacuum,
        gap=gap,
        direction=gb.direction,
        delete_layer=delete_layer,
        tol=tolerance,
        to_primitive=primitive,
    )
    asestruct = AseAtomsAdaptor().get_atoms(structure=gb_struct)
    if graph is not None:
        data = _generate_atomic_sample_data(asestruct, sdict, repeat)
        sample = AtomicScaleSample(**data)

        try:
            gb_inb = GrainBoundary()
            gb_inb.create_grain_boundary(axis=axis, sigma=sigma, gb_plane=gb_plane)
            gb_type = gb_inb.find_gb_character()
        except Exception as e:
            # If GB character detection fails, default to generic GB
            gb_type = None

        if gb_type is None:
            from atomrdf.datamodels.defects.grainboundary import (
                GrainBoundary as GBObject,
            )

            gb_name = "grain_boundary"
        elif gb_type == "Tilt":
            from atomrdf.datamodels.defects.grainboundary import (
                TiltGrainBoundary as GBObject,
            )

            gb_name = "tilt_grain_boundary"
        elif gb_type == "Twist":
            from atomrdf.datamodels.defects.grainboundary import (
                TwistGrainBoundary as GBObject,
            )

            gb_name = "twist_grain_boundary"
        elif gb_type == "Symmetric Tilt":
            from atomrdf.datamodels.defects.grainboundary import (
                SymmetricalTiltGrainBoundary as GBObject,
            )

            gb_name = "symmetric_tilt_grain_boundary"
        elif gb_type == "Mixed":
            from atomrdf.datamodels.defects.grainboundary import (
                MixedGrainBoundary as GBObject,
            )

            gb_name = "mixed_grain_boundary"

        datadict = GBObject.template()
        datadict["sigma"] = sigma  # Use the input parameter
        datadict["rotation_axis"] = axis
        datadict["plane"] = gb_plane
        datadict["misorientation_angle"] = gb.theta[0]  # gb is the AIMSGB object
        setattr(sample, gb_name, GBObject(**datadict))

        sample.to_graph(graph)
        asestruct.info["id"] = sample.id

    return asestruct


def vacancy(
    name,
    indices=None,
    no_of_vacancies=None,
    crystalstructure: str = None,
    a: float = None,
    b: float = None,
    c: float = None,
    *,
    alpha: float = None,
    covera: float = None,
    u: float = None,
    orthorhombic: bool = False,
    cubic: bool = False,
    basis=None,
    repeat: int = 1,
    graph=None,
):

    if indices is None and no_of_vacancies is None:
        raise ValueError("Either indices or no_of_vacancies must be provided")

    if isinstance(name, Atoms):
        # just delet an atom
        if indices is None:
            indices = np.random.choice(range(len(name)), no_of_vacancies, replace=False)
        atoms = name.copy()
        indices = np.sort(indices)
        for index in indices[::-1]:
            del atoms[index]

        # atoms object is provided
        if graph is not None:
            # this means that old system was already linked to a graph
            if "id" in system.info.keys():
                initial_sample_id = system.info["id"]
                # we recreate the sample
                sample = AtomicScaleSample.from_graph(graph, system.info["id"])
                # ok but we have deleted an atom, so we need to update the sample
                sample.update_attributes(atoms)
                # now update the atom attributes
                sample.vacancy = Vacancy(
                    **{
                        "concentration": no_of_vacancies / len(atoms),
                        "number": no_of_vacancies,
                    }
                )
                sample.to_graph(graph)
                final_sample_id = sample.id
                # now we can add the activity
                data = {
                    "initial_sample": initial_sample_id,
                    "final_sample": final_sample_id,
                }
                activity = ops.DeleteAtom(**data)
                activity.to_graph(graph)

            # this means that old system was not linked to a graph
            # the user gave us a graph, so we create a new sample
            else:
                data = _generate_atomic_sample_data(
                    atoms,
                )
                sample = AtomicScaleSample(**data)
                sample.vacancy = Vacancy(
                    **{
                        "concentration": no_of_vacancies / len(atoms),
                        "number": no_of_vacancies,
                    }
                )
                sample.to_graph(graph)

        return atoms

    else:
        atoms, sdict = bulk(
            name,
            crystalstructure=crystalstructure,
            a=a,
            b=b,
            c=c,
            alpha=alpha,
            covera=covera,
            u=u,
            orthorhombic=orthorhombic,
            cubic=cubic,
            basis=basis,
            repeat=repeat,
            get_metadata=True,
        )

        if indices is None:
            indices = np.random.choice(
                range(len(atoms)), no_of_vacancies, replace=False
            )
        indices = np.sort(indices)

        remaining_atoms = len(atoms) - len(indices)
        if remaining_atoms <= 0:
            raise ValueError(
                f"Number of vacancies {len(indices)} is greater than or equal to the number of atoms {len(atoms)}"
            )

        for index in indices[::-1]:
            del atoms[index]

        if graph is not None:
            data = _generate_atomic_sample_data(atoms, sdict, repeat)
            sample = AtomicScaleSample(**data)
            vacancy_data = Vacancy(
                **{
                    "concentration": no_of_vacancies / len(atoms),
                    "number": no_of_vacancies,
                }
            )
            sample.vacancy = vacancy_data
            sample.to_graph(graph)
            atoms.info["id"] = sample.id
        return atoms


def interstitial(
    atoms,
    element,
    void_type="tetrahedral",
    number=1,
    a=None,
    threshold=0.01,
    graph=None,
):
    """
    Create interstitial defects by adding atoms at void positions.

    Parameters
    ----------
    atoms : ase.Atoms or str
        Either an ASE Atoms object or element symbol for creating bulk structure
    element : str or list
        Element symbol(s) for interstitial atom(s)
    void_type : str, optional
        Type of void position: 'tetrahedral' or 'octahedral'. Default is 'tetrahedral'
    number : int, optional
        Number of interstitial atoms to add. Default is 1
    a : float, optional
        Lattice constant. Required for octahedral voids
    threshold : float, optional
        Threshold for finding octahedral positions. Default is 0.01
    graph : KnowledgeGraph, optional
        KnowledgeGraph to add the structure to

    Returns
    -------
    ase.Atoms
        Structure with interstitial defects
    """
    try:
        from pyscal3 import System
    except ImportError:
        raise ImportError("This function requires the pyscal3 package to be installed")

    if not isinstance(atoms, Atoms):
        raise ValueError("atoms must be an ASE Atoms object")

    atoms_copy = atoms.copy()
    sys = System(atoms_copy, format="ase")
    element = np.atleast_1d(element)

    # Ensure we have the right number of elements
    if len(element) < number:
        element = np.repeat(element, number)[:number]
    elif len(element) > number:
        element = element[:number]

    if void_type == "tetrahedral":
        sys.find.neighbors(method="voronoi", cutoff=0.1)
        verts = sys.unique_vertices

        if len(verts) == 0:
            raise ValueError("No tetrahedral void positions found")

        if len(verts) < number:
            raise ValueError(
                f"Only {len(verts)} tetrahedral positions found, but {number} requested"
            )

        randindex = np.random.choice(len(verts), size=number, replace=False)
        randpos = np.array(verts)[randindex]

    elif void_type == "octahedral":
        # Try to get lattice constant from graph if available
        if a is None and graph is not None and "id" in atoms_copy.info:
            try:
                sample_id = atoms_copy.info["id"]
                sample = AtomicScaleSample.from_graph(graph, sample_id)
                if (
                    sample.material is not None
                    and sample.material.crystal_structure is not None
                    and sample.material.crystal_structure.unit_cell is not None
                    and sample.material.crystal_structure.unit_cell.lattice_parameter
                    is not None
                ):
                    a = sample.material.crystal_structure.unit_cell.lattice_parameter[0]
            except:
                pass

        if a is None:
            raise ValueError("Please provide lattice constant a for octahedral voids")

        cutoff = a + threshold * 2
        sys.find.neighbors(method="cutoff", cutoff=cutoff)
        octa_pos = []

        for count, dist in enumerate(sys.atoms.neighbors.distance):
            diffs = np.abs(np.array(dist) - a)
            indices = np.where(diffs < threshold)[0]
            vector = np.array(sys.atoms["diff"][count])[indices]
            vector = sys.atoms.positions[count] + vector / 2
            for vect in vector:
                vect = sys.modify.remap_position_to_box(vect)
                octa_pos.append(vect)

        if len(octa_pos) == 0:
            raise ValueError("No octahedral void positions found")

        octa_pos = np.unique(octa_pos, axis=0)

        if len(octa_pos) < number:
            raise ValueError(
                f"Only {len(octa_pos)} octahedral positions found, but {number} requested"
            )

        randindex = np.random.choice(len(octa_pos), size=number, replace=False)
        randpos = octa_pos[randindex]
    else:
        raise ValueError(
            f'Unknown void type: {void_type}. Use "tetrahedral" or "octahedral"'
        )

    # Add interstitial atoms
    sys = sys.add_atoms({"positions": randpos, "species": element})
    new_atoms = sys.write.ase()

    # Calculate defect properties
    no_of_interstitials = len(randpos)
    original_natoms = len(atoms_copy)
    conc_of_interstitials = no_of_interstitials / original_natoms

    # Update graph if provided
    if graph is not None:
        # Generate new sample data
        data = _generate_atomic_sample_data(new_atoms)
        sample = AtomicScaleSample(**data)

        # Add interstitial defect metadata
        sample.interstitial = Interstitial(
            concentration=conc_of_interstitials,
            number=no_of_interstitials,
        )
        sample.to_graph(graph)
        new_atoms.info["id"] = sample.id

        # If there was a previous sample, add activity to track the operation
        if "id" in atoms_copy.info:
            initial_sample_id = atoms_copy.info["id"]
            final_sample_id = sample.id
            data = {
                "input_sample": initial_sample_id,
                "output_sample": final_sample_id,
            }
            activity = ops.AddAtom(**data)
            activity.to_graph(graph)

    return new_atoms


def substitutional(
    atoms,
    element,
    number=1,
    indices=None,
    graph=None,
):
    """
    Create substitutional defects by replacing atoms with different element(s).

    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object
    element : str
        Element symbol for substitutional atom
    number : int, optional
        Number of substitutions to make. Default is 1
    indices : list or array, optional
        Specific atom indices to substitute. If None, random atoms are chosen
    graph : KnowledgeGraph, optional
        KnowledgeGraph to add the structure to

    Returns
    -------
    ase.Atoms
        Structure with substitutional defects
    """
    if not isinstance(atoms, Atoms):
        raise ValueError("atoms must be an ASE Atoms object")

    atoms_copy = atoms.copy()
    species = atoms_copy.get_chemical_symbols()

    # Select atoms to substitute
    if indices is None:
        if number > len(species):
            raise ValueError(
                f"Number of substitutions {number} exceeds number of atoms {len(species)}"
            )
        indices = np.random.choice(len(species), size=number, replace=False)
    else:
        indices = np.atleast_1d(indices)
        number = len(indices)
        if number > len(species):
            raise ValueError(
                f"Number of indices {number} exceeds number of atoms {len(species)}"
            )

    # Perform substitution
    for idx in indices:
        species[idx] = element
    atoms_copy.set_chemical_symbols(species)

    # Calculate defect properties
    no_of_substitutions = number
    conc_of_substitutions = no_of_substitutions / len(atoms_copy)

    # Update graph if provided
    if graph is not None:
        # Generate new sample data
        data = _generate_atomic_sample_data(atoms_copy)
        sample = AtomicScaleSample(**data)

        # Add substitutional defect metadata
        sample.substitutional = Substitutional(
            concentration=conc_of_substitutions,
            number=no_of_substitutions,
        )
        sample.to_graph(graph)
        atoms_copy.info["id"] = sample.id

        # If there was a previous sample, add activity to track the operation
        if "id" in atoms.info:
            initial_sample_id = atoms.info["id"]
            final_sample_id = sample.id
            data = {
                "input_sample": initial_sample_id,
                "output_sample": final_sample_id,
            }
            activity = ops.SubstituteAtom(**data)
            activity.to_graph(graph)

    return atoms_copy
