# atomRDF

> [!NOTE]
> `atomRDF` was previously called `pyscal-rdf`. 

[![codecov](https://codecov.io/gh/pyscal/atomRDF/graph/badge.svg?token=X7S3TP2MP6)](https://codecov.io/gh/pyscal/atomRDF)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/atomrdf.svg)](https://anaconda.org/conda-forge/atomrdf)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/atomrdf?label=pypi&link=https%3A%2F%2Fpypi.org%2Fproject%2Fatomrdf%2F)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10973374.svg)](https://doi.org/10.5281/zenodo.10973374)

`atomRDF` is a Python tool for ontology-based creation, manipulation, and querying of atomic-scale structures and simulation workflows.

The data model is implemented as [Pydantic](https://docs.pydantic.dev/) classes generated from the
[Conceptual Dictionary for Computational Materials Science](https://github.com/OCDO/conceptual_dictionary)
ontologies maintained by [OCDO](https://github.com/OCDO):

- **CMSO** — Computational Material Sample Ontology
- **CDCO** — Crystallographic Defect Core Ontology
- **PODO**, **PLDO**, **LDO** — point-, planar-, and line-defect ontologies
- **ASMO** — Atomistic Simulation Methods Ontology

This makes every object created by `atomRDF` round-trippable between Python, JSON/YAML and RDF, with
ontology-conformant semantics out of the box.

More details coming soon...


## Acknowledgements
This work is supported by the [NFDI-Matwerk](https://nfdi-matwerk.de/) consortia.

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under the National Research Data Infrastructure – NFDI 38/1 – project number 460247524

