import os

from pyscal_rdf.rdfutils import convert_to_dict
from pyscal.core import System
from pyscal.crystal_structures import Structure
from pyscal_rdf.json_io import write_file
from ase.io import read

def test_poscar():
    sys = System("tests/al_data/Al.poscar", format="poscar")
    convert_to_dict(sys)
    write_file("dump", convert_to_dict(sys))
    assert os.path.exists("dump.json")
    
def test_lammps():
    sys = System("tests/al_data/Al.dump")
    convert_to_dict(sys)
    write_file("dump2", convert_to_dict(sys))
    assert os.path.exists("dump2.json")

def test_cif():
    aseobj = read("tests/al_data/Al.cif", format="cif")
    sys = System(aseobj, format="ase")
    convert_to_dict(sys)
    write_file("dump3", convert_to_dict(sys))
    assert os.path.exists("dump3.json")
    
    
    
    