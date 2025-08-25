import numpy as np
from atomrdf.datamodels.structure import AtomicScaleSample
from atomrdf.build.bulk import _generate_atomic_sample_data


def repeat(system, repetitions, graph=None):
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
    new_system = system.copy()
    new_system = new_system.repeat(repetitions)

    # this means that the system is linked to a graph
    if graph is not None:
        if "id" in system.info.keys():
            # we recreate the sample
            sample = AtomicScaleSample.from_graph(graph, system.info["id"])
            # now update the atom attributes
            sample.update_attributes(new_system, repetitions)
        else:
            data = _generate_atomic_sample_data(new_system, repeat=repetitions)
            sample = AtomicScaleSample(**data)

        # now serialize the sample to the graph
        sample.to_graph(graph)
        # ok this adds a complete new sample - no info about the original, which is ok

    return new_system


def rotate(system, rotation_vectors, graph=None, label=None):
    try:
        import atomman as am
        import atomman.unitconvert as uc
    except ImportError:
        raise ImportError("This function requires the atomman package to be installed")

    system = am.load("ase_Atoms", system)

    # now rotate with atomman
    system = system.rotate(rotation_vectors)

    # get back ase
    new_system = am.dump("ase_Atoms", system)

    if graph is not None:
        if "id" in system.info.keys():
            # we recreate the sample
            sample = AtomicScaleSample.from_graph(graph, system.info["id"])
            # now update the atom attributes
            sample.update_attributes(new_system)
        else:
            data = _generate_atomic_sample_data(
                new_system,
            )
            sample = AtomicScaleSample(**data)

        # now serialize the sample to the graph
        sample.to_graph(graph)

    return new_system


def translate(
    system,
    translation_vector,
    graph=None,
):
    new_system = system.copy()
    new_system.translate(translation_vector)

    if graph is not None:
        if "id" in system.info.keys():
            # we recreate the sample
            sample = AtomicScaleSample.from_graph(graph, system.info["id"])
            # now update the atom attributes
            sample.update_attributes(new_system)
        else:
            data = _generate_atomic_sample_data(
                new_system,
            )
            sample = AtomicScaleSample(**data)

        # now serialize the sample to the graph
        sample.to_graph(graph)

    return new_system


def shear(
    system,
    shear_matrix,
    graph=None,
):
    new_system = system.copy()
    cell = new_system.cell.copy()
    newcell = shear_matrix @ cell
    new_system.set_cell(newcell, scale_atoms=True)

    if graph is not None:
        if "id" in system.info.keys():
            # we recreate the sample
            sample = AtomicScaleSample.from_graph(graph, system.info["id"])
            # now update the atom attributes
            sample.update_attributes(new_system)
        else:
            data = _generate_atomic_sample_data(
                new_system,
            )
            sample = AtomicScaleSample(**data)

        # now serialize the sample to the graph
        sample.to_graph(graph)

    return new_system
