from pyscal_rdf.network.network import OntologyNetwork
import os

#prov = OntologyNetwork('pyscal_rdf/data/prov.rdf', delimiter='#')
def read_ontology():
	cmso = OntologyNetwork(os.path.join(os.path.dirname(__file__), "data/cmso.owl"))
	pldo = OntologyNetwork(os.path.join(os.path.dirname(__file__), "data/pldo.owl"))
	podo = OntologyNetwork(os.path.join(os.path.dirname(__file__), "data/podo.owl"))

	combo = cmso + pldo + podo

	combo.add_path(('cmso:Material', 'cmso:hasDefect', 'pldo:PlanarDefect'))
	combo.add_path(('cmso:Material', 'cmso:hasDefect', 'podo:Vacancy'))
	combo.add_path(('cmso:SimulationCell', 'podo:hasVacancyConcentration', 'float'))
	combo.add_path(('cmso:SimulationCell', 'podo:hasNumberOfVacancies', 'int'))
	combo.add_path(('cmso:SimulationCell', 'prov:wasDerivedFrom', 'cmso:SimulationCell'))
	return combo