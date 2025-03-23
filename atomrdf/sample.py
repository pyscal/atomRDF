"""
Sample

In this module a new Sample class is defined.
"""
from pyscal3.atoms import AttrSetter
from atomrdf.namespace import CMSO, PLDO, PODO, ASMO, PROV, Literal, UNSAFEASMO
from rdflib import RDFS, Namespace, RDF, URIRef
import numpy as np
import uuid
import re
import json

MATH = Namespace("http://purls.helmholtz-metadaten.de/asmo/")

class Sample:
    def __init__(self, name, sample_id, graph):
        self._name = name
        self._sample_id = sample_id
        self._graph = graph

        #create attributes
        self.properties = AttrSetter()
        mapdict = {
            'SimulationCellVolume': self._volume,
            'NumberOfAtoms': self._no_of_atoms,
            'SimulationCellLengths': self._simulation_cell_lengths,
            'SimulationCellLengthX': self._simulation_cell_length_x,
            'SimulationCellLengthY': self._simulation_cell_length_y,
            'SimulationCellLengthZ': self._simulation_cell_length_z,
            'SimulationCellRepetitionX': self._simulation_cell_rep_x,
            'SimulationCellRepetitionY': self._simulation_cell_rep_y,
            'SimulationCellRepetitionZ': self._simulation_cell_rep_z,
            'SimulationCellAngleAlpha': self._simulation_cell_angle_alpha,
            'SimulationCellAngleBeta': self._simulation_cell_angle_beta,
            'SimulationCellAngleGamma': self._simulation_cell_angle_gamma
        }
        self.properties._add_attribute(mapdict)

        labels, props = self._input_properties
        self.inputs = AttrSetter()
        mapdict = {}
        for key, value in zip(labels, props):
            if key is not None:
                mapdict[key] = value
        self.inputs._add_attribute(mapdict)

        labels, props = self._output_properties
        self.outputs = AttrSetter()
        mapdict = {}
        for key, value in zip(labels, props):
            if key is not None:
                mapdict[key] = value
        self.outputs._add_attribute(mapdict)

    def __str__(self):
        return f"{self._name}"
    
    def __repr__(self):
        return f"{self._name}"
    
    @property
    def structure(self):
        return self._graph.get_system_from_sample(self._sample_id)

    @property
    def _volume(self):
        """
        Volume of the sample, PhysicalQuantity
        """
        #get the simulation cell
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        #get the volume, these are volume objects now!
        volume = self._graph.value(simcell, CMSO.hasVolume)
        volume_value = self._graph.value(volume, ASMO.hasValue)

        #find there are already volume attached as calculated property
        inps = [k[2] for k in self._graph.triples((self._sample_id, CMSO.hasCalculatedProperty, None))]
        #initially query if there is property with label volume
        labels = [self._graph.value(inp, RDFS.label) for inp in inps]
        if not "SimulationCellVolume" in labels:
            #add this as a calculated property
            self._graph.add((self._sample_id, CMSO.hasCalculatedProperty, volume))
            parent = volume            
        else:
            parent = inps[labels.index("SimulationCellVolume")]
        return Property(volume_value.toPython(), graph=self._graph, parent=parent, unit='ANGSTROM3', sample_parent=self._sample_id)
    
    @property
    def _no_of_atoms(self):
        no_atoms = self._graph.value(self._sample_id, CMSO.hasNumberOfAtoms)
        inps = [k[2] for k in self._graph.triples((self._sample_id, CMSO.hasCalculatedProperty, None))]
        labels = [self._graph.value(inp, RDFS.label) for inp in inps]
        if not "NumberOfAtoms" in labels:
            parent = self._graph.create_node(URIRef(f'property:{uuid.uuid4()}'), UNSAFEASMO.InputParameter)
            self._graph.add((parent, ASMO.hasValue, Literal(no_atoms.toPython())))
            self._graph.add((parent, RDFS.label, Literal("NumberOfAtoms")))
            self._graph.add((self._sample_id, CMSO.hasCalculatedProperty, parent))
        else:
            parent = inps[labels.index("NumberOfAtoms")]
        return Property(no_atoms.toPython(), graph=self._graph, parent=parent, sample_parent=self._sample_id)

    @property
    def _simulation_cell_lengths(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        simcelllength = self._graph.value(simcell, CMSO.hasLength)
        x = self._graph.value(simcelllength, CMSO.hasLength_x).toPython()
        y = self._graph.value(simcelllength, CMSO.hasLength_y).toPython()
        z = self._graph.value(simcelllength, CMSO.hasLength_z).toPython()
        return Property([x, y, z], graph=self._graph, parent=simcelllength, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_length_x(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        simcelllength = self._graph.value(simcell, CMSO.hasLength)
        x = self._graph.value(simcelllength, CMSO.hasLength_x).toPython()
        return Property(x, graph=self._graph, parent=simcelllength, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_length_y(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        simcelllength = self._graph.value(simcell, CMSO.hasLength)
        y = self._graph.value(simcelllength, CMSO.hasLength_y).toPython()
        return Property(y, graph=self._graph, parent=simcelllength, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_length_z(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        simcelllength = self._graph.value(simcell, CMSO.hasLength)
        z = self._graph.value(simcelllength, CMSO.hasLength_z).toPython()
        return Property(z, graph=self._graph, parent=simcelllength, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_angle_alpha(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        angle = self._graph.value(simcell, CMSO.hasAngle)
        alpha = self._graph.value(angle, CMSO.hasAngle_alpha)
        alpha = alpha.toPython() if alpha is not None else None
        return Property(alpha, graph=self._graph, parent=angle, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_angle_beta(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        angle = self._graph.value(simcell, CMSO.hasAngle)
        beta = self._graph.value(angle, CMSO.hasAngle_beta)
        beta = beta.toPython() if beta is not None else None
        return Property(beta, graph=self._graph, parent=angle, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_angle_gamma(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        angle = self._graph.value(simcell, CMSO.hasAngle)
        gamma = self._graph.value(angle, CMSO.hasAngle_gamma)
        gamma = gamma.toPython() if gamma is not None else None
        return Property(gamma, graph=self._graph, parent=angle, sample_parent=self._sample_id)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

    @property
    def _simulation_cell_rep_x(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        x = self._graph.value(simcell, CMSO.hasRepetition_x)
        x = x.toPython() if x is not None else 1
        return Property(x, graph=self._graph, parent=simcell, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_rep_y(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        y = self._graph.value(simcell, CMSO.hasRepetition_y)
        y = y.toPython() if y is not None else 1
        return Property(y, graph=self._graph, parent=simcell, sample_parent=self._sample_id)
    
    @property
    def _simulation_cell_rep_z(self):
        simcell = self._graph.value(self._sample_id, CMSO.hasSimulationCell)
        z = self._graph.value(simcell, CMSO.hasRepetition_z)
        z = z.toPython() if z is not None else 1
        return Property(z, graph=self._graph, parent=simcell, sample_parent=self._sample_id)

    @property
    def _input_properties(self):
        activity = self._graph.value(self._sample_id, PROV.wasGeneratedBy)
        if activity is not None:
            inps = [k[2] for k in self._graph.triples((activity, ASMO.hasInputParameter, None))]
            labels = [self._graph.value(inp, RDFS.label) for inp in inps]
            labels = [label if label is None else label.toPython() for label in labels]
            values = [self._graph.value(inp, ASMO.hasValue) for inp in inps]
            units = [self._graph.value(inp, ASMO.hasUnit) for inp in inps]
            units = [unit if unit is None else unit.toPython().split('/')[-1] for unit in units]
            props = []
            for count, value in enumerate(values):
                props.append(Property(value.toPython(), unit=units[count], graph=self._graph, parent=inps[count], sample_parent=self._sample_id))
            return labels, props
        return [], []
    
    @property
    def _output_properties(self):
        inps = [k[2] for k in self._graph.triples((self._sample_id, CMSO.hasCalculatedProperty, None))]
        labels = [self._graph.value(inp, RDFS.label) for inp in inps]
        labels = [label if label is None else label.toPython() for label in labels]
        values = [self._graph.value(inp, ASMO.hasValue) for inp in inps]
        units = [self._graph.value(inp, ASMO.hasUnit) for inp in inps]
        units = [unit if unit is None else unit.toPython().split('/')[-1] for unit in units]
        props = []
        for count, value in enumerate(values):
            props.append(Property(value.toPython(), unit=units[count], graph=self._graph, parent=inps[count], sample_parent=self._sample_id))
        return labels, props

class Property:
    def __init__(self, value, unit=None, graph=None, parent=None, sample_parent=None):
        self._value = self._clean_value(value)
        self._unit = unit
        self._graph = graph
        self._parent = parent
        self._label = None
        self._sample_parent = sample_parent
    
    def __getitem__(self, index):
        if isinstance(self._value, (list, np.ndarray)):
            if isinstance(index, slice):
                value = self._value[index.start:index.stop:index.step]
            else:
                value = self._value[index]
            return Property(value, unit=self._unit, graph=self._graph, parent=self._parent, sample_parent=self._sample_parent)
            
        else:
            raise TypeError('This property is not a list or array, therefore not subscriptable.')
    

    def _clean_value(self, value):
        if isinstance(value, str):
            if (value[0] == '[') and (value[-1] == ']'):
                value = np.array(json.loads(value))
        if isinstance(value, (list, np.ndarray)):
            value = np.array(value)
        return value
    
    def __repr__(self):
        if self._unit is not None:
            return f"{self._value} {self._unit}"
        return f"{self._value}"
    
    #
    @property
    def value(self):
        return self._value
    
    @property
    def label(self):
        if self._graph is not None:
            label = self._graph.value(self._parent, RDFS.label)
            if label is not None:
                return label.toPython()
        return self._label

    @label.setter
    def label(self, value):
        self._label = value
        if self._graph is not None:
            if self._parent is not None:
                self._graph.remove((self._parent, RDFS.label, None))
                self._graph.add((self._parent, RDFS.label, Literal(value)))
    
    def _create_label(self, v1, v2, operation):
        if isinstance(v1, Property):
            v1 = v1.label
        else:
            v1 = f'v{str(v1)}'
        if isinstance(v2, Property):
            v2 = v2.label
        else:
            v2 = f'v{str(v2)}'
        return f'({v1}{operation}{v2})'

    def _declass(self, item):
        if isinstance(item, Property):
            return item.value
        else:
            return item

    def _wrap(self, value):
        if isinstance(value, Property):
            return value._parent
        else:
            return Literal(value)

    #create node with units
    def _create_node(self, res):
        parent = self._graph.create_node(URIRef(f'property:{uuid.uuid4()}'), CMSO.CalculatedProperty)
        self._graph.add((parent, ASMO.hasValue, Literal(res)))
        if self._unit is not None:
            self._graph.add((parent, ASMO.hasUnit, URIRef(f"http://qudt.org/vocab/unit/{self._unit}")))
        return parent

    def update_physical_quantity(self, ontoterm):
        if self._graph is not None:
            if self._parent is not None:
                self._graph.remove((self._parent, RDF.type, None))
            #remove the old type
            self._graph.add((self._parent, RDF.type, ontoterm.URIRef))

    #overloaded operations
    def __add__(self, value):
        res = self._value + self._declass(value)
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent) 
        res_prop.label = self._create_label(self, value, '+')
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Addition))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasAddend, self._wrap(value)))
            self._graph.add((operation, MATH.hasAddend, self._wrap(self)))
            self._graph.add((operation, MATH.hasSum, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop
    
    def __radd__(self, value):
        return self.__add__(value)
    
    def __sub__(self, value):
        res = self._value - self._declass(value)
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '-') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Subtraction))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasMinuend, self._wrap(self)))
            self._graph.add((operation, MATH.hasSubtrahend, self._wrap(value)))
            self._graph.add((operation, MATH.hasDifference, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop

    def __rsub__(self, value):
        res =  self._declass(value) - self._value
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '-') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Subtraction))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasMinuend, self._wrap(value)))
            self._graph.add((operation, MATH.hasSubtrahend, self._wrap(self)))
            self._graph.add((operation, MATH.hasDifference, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop  
    
    def __mul__(self, value):
        res = self._value * self._declass(value)
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '*') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Multiplication))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasFactor, self._wrap(self)))
            self._graph.add((operation, MATH.hasFactor, self._wrap(value)))
            self._graph.add((operation, MATH.hasProduct, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop
    
    def __rmul__(self, value):
        return self.__mul__(value)

    def __truediv__(self, value):
        res = self._value / self._declass(value)
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '/') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Division))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasDivisor, self._wrap(self)))
            self._graph.add((operation, MATH.hasDividend, self._wrap(value)))
            self._graph.add((operation, MATH.hasQuotient, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop

    def __pow__(self, value):
        res = self._value ** self._declass(value)
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '^') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Exponentiation))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasBase, self._wrap(self)))
            self._graph.add((operation, MATH.hasExponent, self._wrap(value)))
            self._graph.add((operation, MATH.hasPower, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop

    def __rpow__(self, value):
        res =  self._declass(value) ** self._value
        parent = self._create_node(res)
        res_prop = Property(res, unit=self._unit, graph=self._graph, parent=parent)
        res_prop.label = self._create_label(self, value, '^') 
        res_prop._sample_parent = self._sample_parent
        if self._graph is not None:
            operation = URIRef(f'operation:{uuid.uuid4()}')
            self._graph.add((operation, RDF.type, MATH.Exponentiation))
            self._graph.add((operation, RDF.type, PROV.Activity))
            self._graph.add((operation, MATH.hasBase, self._wrap(value)))
            self._graph.add((operation, MATH.hasExponent, self._wrap(self)))
            self._graph.add((operation, MATH.hasPower, self._wrap(res_prop)))
            self._graph.add((parent, MATH.wasCalculatedBy, operation))
        return res_prop
    
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
        return Property(-self._value, unit=self._unit, graph=self._graph, parent=self._parent)
    
    def __abs__(self):
        return Property(abs(self._value), unit=self._unit, graph=self._graph, parent=self._parent)
    
    def __round__(self, n):
        return Property(round(self._value, n), unit=self._unit, graph=self._graph, parent=self._parent)