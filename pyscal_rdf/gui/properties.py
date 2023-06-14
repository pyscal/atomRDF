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
    "GB sigma value ": {"onto": "Sigma", "dtype": "integer"},
    "GB Miller indices ": {"onto": "MillerIndices", "dtype": "string"},
    "GB rotation axis x ": {"onto": "RotationAxis_x", "dtype": "float"},
    "GB rotation axis y ": {"onto": "RotationAxis_y", "dtype": "float"},
    "GB rotation axis z ": {"onto": "RotationAxis_z", "dtype": "float"},
    "GB misorientation angle ": {"onto": "Angle", "dtype": "float"},
    "vacancy concentration ": {"onto": "VacancyConcentration", "dtype": "float"},
    "number of vacancy ": {"onto": "NumberOfVacancy", "dtype": "integer"},
}

classes = {
    "simulation cell": {"onto": "SimulationCell"},
}

options = list(dataprops.keys())
#options = list(dataprops.keys()) + list(classes.keys())