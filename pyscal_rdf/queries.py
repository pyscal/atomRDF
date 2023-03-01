from rdflib import URIRef, Namespace, BNode, Literal, XSD
from rdflib.namespace import RDF

CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")

class Query:
    def __init__(self):
        self.sparql = Sparql()
        self.python = Python()

class Sparql:
    def __init__(self):
        pass
    
    def sample_by_altname(self, g, altname):
        query = """
        PREFIX cmso: <https://purls.helmholtz-metadaten.de/cmso/>
        SELECT DISTINCT ?sample
        WHERE {
            ?crystalstructure cmso:hasAltName ?altname .
            ?material cmso:hasStructure ?crystalstructure .
            ?sample cmso:hasMaterial ?material .

        FILTER (?altname="%s"^^xsd:string)
        }"""%(altname)
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
                
            
        

    
