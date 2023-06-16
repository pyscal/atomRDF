import pytest
from pyscal.core import System
from pyscal_rdf.rdfutils import convert_to_dict
import pyscal_rdf.properties as prp

def test_convert_to_dict():
    sys = System.from_structure("bcc", repetitions=[1,1,1], element="Fe")
    info = convert_to_dict(sys)
    assert(info["ChemicalCompositionElement"] == list(sys.composition.keys()))
    assert(info["ChemicalCompositionRatio"] == [val for key, val in sys.composition.items()])    
    assert(info["CellVolume"] == sys.volume)
    assert(info["NumberOfAtoms"] == sys.natoms)
    
    box_dims = sys.box_dimensions
    assert(info["SimulationCellLengthX"] == box_dims[0])
    assert(info["SimulationCellLengthY"] == box_dims[1])
    assert(info["SimulationCellLengthZ"] == box_dims[2])

    assert(info["SimulationCellVectorA"] == sys.box[0])
    assert(info["SimulationCellVectorB"] == sys.box[1])
    assert(info["SimulationCellVectorC"] == sys.box[2])
    
    assert(info["SimulationCellAngleAlpha"] == prp.get_angle(sys.box[0], sys.box[1]))
    assert(info["SimulationCellAngleBeta"] == prp.get_angle(sys.box[1], sys.box[2]))
    assert(info["SimulationCellAngleGamma"] == prp.get_angle(sys.box[2], sys.box[0]))    
    
    #Atom properties
    #--------------------------------------------------
    if sys.atoms.species[0] is not None:
        assert(info["Element"] == sys.atoms.species)
    else:
        assert(info["Element"] == sys.atoms.types)
    assert(info["Coordination"] == prp.get_coordination(sys))
    assert(info["Positions"] == sys.atoms.positions)
    
    #unit cell things
    #--------------------------------------------------
    #Unit cell properties
    assert(info["LatticeParameter"] == sys.atoms._lattice_constant)
    
    #do only for unit cell
    if sys._structure_dict is None:
        pass
        #maybe try to find space group with warnings!
    else:
        symbol, number = prp.get_space_group(sys)
        assert(info["SpaceGroupSymbol"] == symbol)
        assert(info["SpaceGroupNumber"] == number)
        assert(info["CrystalStructureName"] == sys.atoms._lattice)
        assert(info["BravaisLattice"] == prp.get_bravais_lattice(sys))
        assert(info["BasisPositions"] == sys._structure_dict['positions'])
        assert(info["BasisOccupancy"] == prp.get_basis(sys))
        assert(info["LatticeVectors"] == prp.get_lattice_vector(sys))
