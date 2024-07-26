from atomrdf import KnowledgeGraph, System
from ase.io import read
import os

def test_qe_write():
    kg = KnowledgeGraph()
    s = System.create.element.Fe(graph = kg)
    s.write.file('tests/pw.si.scf.in',         
            format='quantum-espresso',
            copy_from='tests/qe_data/pw.si.scf.in',
            pseudo_files=['tests/qe_data/Si.pbe-n-rrkjus_psl.1.0.0.UPF'])
    assert os.path.exists('tests/pw.si.scf.in')
    struct = read('tests/pw.si.scf.in', format='espresso-in')
    assert len(struct) == 2

def test_qe_read():
    kg = KnowledgeGraph()
    job = ('tests/qe_data/pw.si.scf.in', 'tests/qe_data/pw.si.scf.out')
    kg.add_workflow(job, workflow_environment='qe')
    assert kg.n_samples == 2