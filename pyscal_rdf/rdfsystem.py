import numpy as np
import pyscal3.core as pc

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

	def __delitem__(self, val):
        if isinstance(val, int):
            val = [val]
        self.delete(indices=list(val))

