import numpy as np
import pyscal3.core as pc
from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PLDO = Namespace("https://purls.helmholtz-metadaten.de/pldo/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")

class System(pc.System):
	def __init__(self, filename = None, 
			format = "lammps-dump", 
        	compressed = False, 
        	customkeys = None):
		super().__init__(filename = filename, 
			format = format, 
        	compressed = compressed, 
        	customkeys = customkeys)
		#this is the sample which will be stored
		self.sample = None
		#the graph object should also be attached
		#for post-processing of structures
		self.graph = None
		self._atom_ids = None

	def __delitem__(self, val):
        if isinstance(val, int):
            val = [val]
        self.delete(indices=list(val))
        #now the graph has to be updated accordingly
        if self.graph is not None:
        	#first annotate graph
        	c = (len(val)/self.natoms)
        	self.graph.add_vacancy(c, number=len(val))
        	#now we need to re-add atoms, so at to remove
        	#deleted ones from the vacancy
        	atoms = list([s[2] for s in self.graph.triples((self.sample, CMSO.hasAtom, None))])
        	#this is the list of atoms in this sample
        	for atom in atoms:
        		self.graph.remove(())

        	self.graph.remove((None, CMSO.hasAtom, None))
        	self.graph.remove((None, None, CMSO.Atom))
        	self.graph.remove((None, CMSO.hasPositionVector, None))
        	self.graph.remove(())
        	self.graph.remove(())
        	self.graph.remove(())
        	self.graph.remove(())
        	self.graph.remove(())

