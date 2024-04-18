import pytest
from atomrdf.encoder import NumpyArrayEncoder
from atomrdf.json_io import write_file
import json
import numpy as np
import os

def test_encoder():
    arr = np.linspace(0, 10, 100)
    data = {"array": arr}
    dumpdata = json.dumps(data, cls=NumpyArrayEncoder)
    assert(np.abs(float(dumpdata.split(',')[-2]) - arr[-2]) < 1E-5)

    enc = NumpyArrayEncoder()
    assert type(enc.default(np.float64(23))) == float
    assert type(enc.default(np.int64(23))) == int
    assert type(enc.default(np.array([23]))) == list

def test_writer(tmp_path):
    arr = np.linspace(0, 10, 100)
    data = {"array": arr}
    write_file("test_json", data)
    assert(os.path.exists("test_json.json"))

