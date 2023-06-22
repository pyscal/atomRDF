from rdflib import URIRef, Namespace, BNode, Literal, XSD
from rdflib.namespace import RDF
import numpy as np

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")

class Query:
    def __init__(self):
        self.sparql = Sparql()
        self.python = Python()

class Sparql:
    def __init__(self):
        pass
    
    def sample_by_latticesystem(self, g, latticesystem):
        query = """
        PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
        PREFIX wiki: <https://www.wikidata.org/wiki/>
        SELECT DISTINCT ?sample
        WHERE {
            ?unitcell cmso:hasBravaisLattice wiki:%s .
            ?structure cmso:hasUnitCell ?unitcell .
            ?material cmso:hasStructure ?structure .
            ?sample cmso:hasMaterial ?material .
        }"""%(latticesystem)
        qres = g.graph.query(query)
        samples = [j[0] for j in qres]
        return samples
    
    def sample_by_sigma(self, g, sigma):
        if isinstance(sigma, int):
            query="""
            PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
            PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            SELECT DISTINCT ?sample
            WHERE {
                ?sample cmso:hasMaterial ?material .
                ?material cmso:hasDefect ?defect .
                ?defect pldo:hasSigmaValue ?sigma .
            FILTER (?sigma="%d"^^xsd:integer)
            }"""%(sigma)
        elif isinstance(sigma, list):
            if not len(sigma) == 2:
                raise ValueError("Range queries support only arrays of length 2")
            sigma = np.sort(np.abs(sigma))
            query="""
            PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
            PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            SELECT DISTINCT ?sample
            WHERE {
                ?sample cmso:hasMaterial ?material .
                ?material cmso:hasDefect ?defect .
                ?defect pldo:hasSigmaValue ?sigma .
            FILTER (?sigma >= "%d"^^xsd:integer && ?sigma <= "%d"^^xsd:integer)
            }"""%(sigma[0], sigma[1])
        qres = g.graph.query(query)
        samples = [j[0] for j in qres]
        return samples          
            
    def sample_by_defect(self, g, defect):
        defect = defect.lower()
        defectdict = {"tilt": "TiltBoundary", 
                      "twist": "TwistBoundary", 
                      "symmetric tilt": "SymmetricTiltBoundary", 
                      "mixed": "MixedBoundary", 
                      "grain boundary": None}
        
        if defect not in defectdict.keys():
            raise ValueError(f'Please choose defect from {defectlist.keys()}')
        if defect != "grain boundary":
            query = """
            PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
            PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            SELECT DISTINCT ?sample
            WHERE {
                ?sample cmso:hasMaterial ?material .
                ?material cmso:hasDefect ?defect .
                ?defect rdf:type pldo:%s .
            }"""%(defectdict[defect])
        else:
            querycomm = ['{?defect rdf:type pldo:%s}'%val for key, val in defectdict.items()]
            querycomm = " UNION ".join(querycomm)
            query = """
            PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
            PREFIX pldo: <https://purls.helmholtz-metadaten.de/pldo/>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            SELECT DISTINCT ?sample
            WHERE {
                ?sample cmso:hasMaterial ?material .
                ?material cmso:hasDefect ?defect .
                %s .
            }"""%(querycomm)
        qres = g.graph.query(query)
        samples = [j[0] for j in qres]
        return samples        
class Python:
    def __init__(self):
        pass
    
    def sample_by_altname(self, g, altname):
        samples = []
        for s in g.graph.triples((None, CMSO.hasAltName, Literal(altname, datatype=XSD.string))):
            for k in g.graph.triples((None, CMSO.hasStructure, s[0])):
                for j in g.graph.triples((None, CMSO.hasMaterial, k[0])):
                    samples.append(j[0])
        return samples
                
            
        

    
