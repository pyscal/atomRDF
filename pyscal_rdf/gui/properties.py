dataprops = {
    "volume": {"onto": "Volume", "dtype": "float"},
    "number of atoms": {"onto": "NumberOfAtoms", "dtype": "integer"},
    "simulation cell length in x ": {"onto": "SimulationCellLength_x", "dtype": "float"},
    "simulation cell length in y ": {"onto": "SimulationCellLength_y", "dtype": "float"},
    "simulation cell length in z ": {"onto": "SimulationCellLength_z", "dtype": "float"},
    "crystal structure alternative name ": {"onto": "CrystalStructureAltName", "dtype": "string"},
    "space group symbol ": {"onto": "SpaceGroupSymbol", "dtype": "string"},
    "space group number ": {"onto": "SpaceGroupNumber", "dtype": "integer"},
    "Bravais lattice ": {"onto": "LatticeSystem", "dtype": "string"},
    "lattice parameter in x ": {"onto": "LatticeParameter_x", "dtype": "float"},
    "lattice parameter in y ": {"onto": "LatticeParameter_y", "dtype": "float"},
    "lattice parameter in z ": {"onto": "LatticeParameter_z", "dtype": "float"},
}

classes = {
    "simulation cell": {"onto": "SimulationCell"},
}

options = list(dataprops.keys())
#options = list(dataprops.keys()) + list(classes.keys())