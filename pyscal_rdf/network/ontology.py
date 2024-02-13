from pyscal_rdf.network.network import OntologyNetwork
import os

#prov = OntologyNetwork('pyscal_rdf/data/prov.rdf', delimiter='#')
def read_ontology():
	#read in ontologies
	file_location = os.path.dirname(__file__).split('/')
	file_location = "/".join(file_location[:-1])

	cmso = OntologyNetwork(os.path.join(file_location,  'data/cmso.owl'))
	pldo = OntologyNetwork(os.path.join(file_location,  'data/pldo.owl'))
	podo = OntologyNetwork(os.path.join(file_location,  'data/podo.owl'))
	asmo = OntologyNetwork(os.path.join(file_location,  'data/asmo.owl'))
	
	#combine them	
	combo = cmso + pldo + podo + asmo

	#add namespaces
	combo.add_namespace('prov', 'http://www.w3.org/ns/prov#')
	combo.add_namespace('rdf', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#')
	combo.add_namespace('rdfs', 'http://www.w3.org/2000/01/rdf-schema#')
	
	#add extra terms for quering
	combo.add_term('http://www.w3.org/ns/prov#Entity', 'class', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#Activity', 'class', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#SoftwareAgent', 'class', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#wasDerivedFrom', 'object_property', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#wasGeneratedBy', 'object_property', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#wasAssociatedWith', 'object_property', delimiter='#')
	combo.add_term('http://www.w3.org/ns/prov#actedOnBehalfOf', 'object_property', delimiter='#')
	combo.add_term('http://www.w3.org/2000/01/rdf-schema#label', 'data_property', delimiter='#', namespace='rdfs')
	combo.add_term('http://www.w3.org/1999/02/22-rdf-syntax-ns#type', 'object_property', delimiter='#', namespace='rdf')
	
	#add paths

	#General fixes
	combo.add_path(('cmso:CrystalStructure', 'cmso:hasAltName', 'string'))
	combo.add_path(('cmso:ChemicalElement', 'cmso:hasSymbol', 'string'))

	#interontology paths
	combo.add_path(('cmso:Material', 'cmso:hasDefect', 'pldo:PlanarDefect'))
	combo.add_path(('cmso:Material', 'cmso:hasDefect', 'podo:Vacancy'))
	combo.add_path(('cmso:SimulationCell', 'podo:hasVacancyConcentration', 'float'))
	combo.add_path(('cmso:SimulationCell', 'podo:hasNumberOfVacancies', 'int'))
	combo.add_path(('cmso:ComputationalSample', 'prov:wasDerivedFrom', 'cmso:ComputationalSample'))
	combo.add_path(('cmso:ComputationalSample', 'rdf:type', 'prov:Entity'))
	combo.add_path(('asmo:StructureOptimization', 'rdf:type', 'prov:Activity'))	
	combo.add_path(('asmo:StructureOptimization', 'prov:wasAssociatedWith', 'prov:SoftwareAgent'))
	combo.add_path(('cmso:ComputationalSample', 'prov:wasGeneratedBy', 'asmo:StructureOptimization'))
	
	#return
	return combo