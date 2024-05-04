"""
Sample

In this module a new Sample class is defined.
"""
from pyscal3.atoms import AttrSetter
from atomrdf.namespace import CMSO, PLDO, PODO

class Sample:
    def __init__(self, name, sample_id, graph):
        self._name = name
        self._sample_id = sample_id
        self._graph = graph

    def __str__(self):
        return f"{self._name}"
    
    def __repr__(self):
        return f"{self._name}"
    
    @property
    def _volume(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        volume = self._graph.value(simcell, CMSO.hasVolume)
        return Property(volume.toPython())
    

class Property:
    def __init__(self, value):
        self._value = value
    
    def __repr__(self):
        return f"{self._value}"
    
    @property
    def value(self):
        return self._value
    
    #overloaded operations
    def __add__(self, other):
        return Property(self._value + other.value)
    
    def __sub__(self, other):
        return Property(self._value - other.value)
    
    def __mul__(self, other):
        return Property(self._value * other.value)

    def __truediv__(self, other):
        return Property(self._value / other.value)
    
    def __eq__(self, other):
        return self._value == other.value

    def __ne__(self, other):
        return self._value != other.value
    
    def __lt__(self, other):
        return self._value < other.value
    
    def __le__(self, other):
        return self._value <= other.value
    
    def __gt__(self, other):
        return self._value > other.value
    
    def __ge__(self, other):
        return self._value >= other.value
    
    def __neg__(self):
        return Property(-self._value)
    
    def __abs__(self):
        return Property(abs(self._value))
    
    def __round__(self, n):
        return Property(round(self._value, n))