import pyscal3.operations.operations as operations

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

