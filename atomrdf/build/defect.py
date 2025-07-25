from platform import system
from atomrdf.build.bulk import bulk, _generate_atomic_sample_data
import numpy as np
from atomrdf.datamodels.defects.dislocation import Dislocation
from atomrdf.datamodels.defects.stackingfault import StackingFault
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.build.buildutils import _declass


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
        datadict["plabe"]["value"] = slip_plane
        datadict["displacement"]["value"] = displ
        setattr(sample, "stacking_fault", StackingFault(**datadict))
        sample.to_graph(graph)

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
        disl_dict["line_direction"]["value"] = dislocation_line
        disl_dict["burgers_vector"]["value"] = burgers_vector
        disl_dict["slip_system"]["slip_direction"]["value"] = slip_direction
        disl_dict["slip_system"]["slip_plane"]["normal"]["value"] = slip_plane
        if disl_name == "mixed_dislocation":
            disl_dict["character_angle"]["value"] = angle_deg

        setattr(sample, disl_name, DislocationObject(**disl_dict))
        sample.to_graph(graph)

    if return_atomman_dislocation:
        return aseatoms, disc
    return aseatoms
