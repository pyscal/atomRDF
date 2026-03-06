import json
import numpy as np
from rdflib import Namespace, URIRef
from rdflib import Literal as RDFLibLiteral
from rdflib import Namespace as RDFLibNamespace


def Literal(value, datatype=None):
    if datatype is not None:
        return RDFLibLiteral(value, datatype=datatype)
    elif isinstance(value, list):
        return RDFLibLiteral(json.dumps(value))
    elif isinstance(value, np.ndarray):
        return RDFLibLiteral(json.dumps(value.tolist()))
    else:
        return RDFLibLiteral(value)


CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
ASMO = Namespace("http://purls.helmholtz-metadaten.de/asmo/")
LDO  = Namespace("http://purls.helmholtz-metadaten.de/cdos/ldo/")
PLDO = Namespace("http://purls.helmholtz-metadaten.de/cdos/pldo/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/cdos/podo/")
CDCO = Namespace("http://purls.helmholtz-metadaten.de/cdos/cdco/")

PROV = Namespace("http://www.w3.org/ns/prov#")
MDO  = Namespace("https://w3id.org/mdo/calculation/")
DCAT = Namespace("http://www.w3.org/ns/dcat#")
