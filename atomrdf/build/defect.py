from atomrdf.build.bulk import bulk


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
    repeat=None,
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

    if structure is not None:
        # create a structure with the info
        input_structure = bulk(
            element,
            structure=structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=noise,
            primitive=primitive,
        )
    elif element is not None:
        if element in element_dict.keys():
            structure = element_dict[element]["structure"]
            lattice_constant = element_dict[element]["lattice_constant"]
        else:
            raise ValueError("Please provide structure")
        input_structure = bulk(
            element,
            structure=structure,
            lattice_constant=_declass(lattice_constant),
            repetitions=repetitions,
            ca_ratio=_declass(ca_ratio),
            noise=noise,
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

    # find dislocation character
    angle = np.dot(dislocation_line, burgers_vector) / (
        np.linalg.norm(dislocation_line) * np.linalg.norm(burgers_vector)
    )
    angle_rad = np.arccos(angle)
    angle_deg = np.degrees(angle_rad)

    disl_dict = {
        "BurgersVector": burgers_vector,
        "SlipPlane": slip_plane,
        "SlipDirection": slip_direction,
        "DislocationLine": dislocation_line,
        "DislocationCharacter": angle_deg,
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
    output_structure.add_property_mappings(
        lattice_constant, mapping_quantity="lattice_constant"
    )
    output_structure.add_property_mappings(
        ca_ratio, mapping_quantity="lattice_constant"
    )

    if return_atomman_dislocation:
        return output_structure, disc
    return output_structure
