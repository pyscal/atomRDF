"""
Sample

In this module a new Sample class is defined.
"""
from pyscal3.atoms import AttrSetter
from atomrdf.namespace import CMSO, PLDO, PODO, ASMO, PROV
from rdflib import RDFS
import numpy as np

class Sample:
    def __init__(self, name, sample_id, graph):
        self._name = name
        self._sample_id = sample_id
        self._graph = graph

        #create attributes
        self.properties = AttrSetter()
        mapdict = {
            'volume': self._volume,
            'no_of_atoms': self._no_of_atoms,
        }
        self.properties._add_attribute(mapdict)

        self.inputs = AttrSetter()
        mapdict = {}
        for key, value in zip(*self._input_properties):
            mapdict[key] = value
        self.inputs._add_attribute(mapdict)

        self.outputs = AttrSetter()
        mapdict = {}
        for key, value in zip(*self._output_properties):
            mapdict[key] = value
        self.outputs._add_attribute(mapdict)
        
    def __str__(self):
        return f"{self._name}"
    
    def __repr__(self):
        return f"{self._name}"
    
    @property
    def _volume(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        volume = self._graph.value(simcell, CMSO.hasVolume)
        return Property(volume.toPython(), graph=self._graph)
    
    @property
    def _no_of_atoms(self):
        no_atoms = self._graph.value(self._sample_id, CMSO.hasNumberOfAtoms)
        return Property(no_atoms.toPython(), graph=self._graph)
    
    @property
    def _input_properties(self):
        activity = self._graph.value(self._sample_id, PROV.wasGeneratedBy)
        if activity is not None:
            inps = [k[2] for k in self._graph.triples((activity, ASMO.hasInputParameter, None))]
            labels = [self._graph.value(inp, RDFS.label) for inp in inps]
            values = [self._graph.value(inp, ASMO.hasValue) for inp in inps]
            units = [self._graph.value(inp, ASMO.hasUnit) for inp in inps]
            units = [unit if unit is None else unit.toPython().split('/')[-1] for unit in units]

            props = [Property(value.toPython(), unit=unit, graph=self._graph) for value, unit in zip(values, units)]
            return props, labels

    @property
    def _output_properties(self):
        inps = [k[2] for k in self._graph.triples((self._sample_id, CMSO.hasCalculatedProperty, None))]
        labels = [self._graph.value(inp, RDFS.label) for inp in inps]
        values = [self._graph.value(inp, ASMO.hasValue) for inp in inps]
        units = [self._graph.value(inp, ASMO.hasUnit) for inp in inps]
        units = [unit if unit is None else unit.toPython().split('/')[-1] for unit in units]

        props = [Property(value.toPython(), unit=unit, graph=self._graph) for value, unit in zip(values, units)]
        return props, labels   

class Property:
    def __init__(self, value, unit=None, graph=None):
        self._value = value
        self._unit = unit
        self._graph = graph
    
    def __repr__(self):
        if self._unit is not None:
            return f"{self._value} {self._unit}"
        return f"{self._value}"
    
    @property
    def value(self):
        return self._value
    
    def _declass(self, item):
        if isinstance(item, Property):
            return item.value
        else:
            return item
        
    #overloaded operations
    def __add__(self, value):
        return Property(self._value + self._declass(value), unit=self._unit, graph=self._graph)
    
    def __sub__(self, value):
        return Property(self._value - self._declass(value), unit=self._unit, graph=self._graph)
    
    def __mul__(self, value):
        return Property(self._value * self._declass(value), unit=self._unit, graph=self._graph)

    def __truediv__(self, value):
        return Property(self._value / self._declass(value), unit=self._unit, graph=self._graph)
    
    def __eq__(self, value):
        return self._value == self._declass(value)

    def __ne__(self, value):
        return self._value != self._declass(value)
    
    def __lt__(self, value):
        return self._value < self._declass(value)
    
    def __le__(self, value):
        return self._value <= self._declass(value)
    
    def __gt__(self, value):
        return self._value > self._declass(value)
    
    def __ge__(self, value):
        return self._value >= self._declass(value)
    
    def __neg__(self):
        return Property(-self._value, unit=self._unit, graph=self._graph)
    
    def __abs__(self):
        return Property(abs(self._value), unit=self._unit, graph=self._graph)
    
    def __round__(self, n):
        return Property(round(self._value, n), unit=self._unit, graph=self._graph)