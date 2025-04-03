from atomrdf.structure import System

def read(
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
        type 2 will be matc_read_structurehed to 'Al'

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
        specify the number _read_structureof repetitions of the unit cell in each direction. Default is None.

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