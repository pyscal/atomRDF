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

## Installation

```bash
pip install atomrdf
# or via conda-forge
conda install -c conda-forge atomrdf
```

Optional features ship as extras &mdash; install only what you need:

```bash
pip install "atomrdf[oxigraph]"           # Oxigraph triple-store backend
pip install "atomrdf[sqlalchemy]"         # SQLAlchemy-backed store
pip install "atomrdf[materials_project]"  # Materials Project lookups (mp-api)
pip install "atomrdf[grainboundary]"      # aimsgb + pymatgen for grain boundaries
pip install "atomrdf[dislocation]"        # atomman for dislocation builders
```

## Quickstart

```python
from atomrdf import KnowledgeGraph
import atomrdf.build as build

# 1. open a knowledge graph (in-memory by default)
kg = KnowledgeGraph()

# 2. build a sample; it is automatically annotated and added to the graph
fe = build.bulk("Fe", cubic=True, graph=kg)

# 3. ask SPARQL questions
results = kg.query("""
    PREFIX cmso: <http://purls.helmholtz-metadaten.de/cmso/>
    SELECT ?sample ?n
    WHERE { ?sample cmso:hasNumberOfAtoms ?n }
""")
print(results)

# 4. persist
kg.write("fe.ttl", format="turtle")
```

See [`examples/`](examples/) for end-to-end notebooks (getting started, grain
boundaries, working with data, defects, SPARQL queries, …) and the full
documentation at <https://atomrdf.pyscal.org>.

## Citing atomRDF

If you use atomRDF in academic work, please cite:

> Guzmán, A. A., Menon, S., Hickel, T., & Sandfeld, S. (2026).
> *Ontology-based knowledge graph infrastructure for interoperable atomistic
> simulation data.* arXiv:2604.06230. <https://arxiv.org/abs/2604.06230>

BibTeX:

```bibtex
@misc{guzman2026ontologybasedknowledgegraphinfrastructure,
  title         = {Ontology-based knowledge graph infrastructure for interoperable atomistic simulation data},
  author        = {Abril Azocar Guzman and Sarath Menon and Tilmann Hickel and Stefan Sandfeld},
  year          = {2026},
  eprint        = {2604.06230},
  archivePrefix = {arXiv},
  primaryClass  = {cs.DB},
  url           = {https://arxiv.org/abs/2604.06230},
}
```

A software citation (Zenodo DOI) is also available via the badge above and
[CITATION.cff](CITATION.cff).


## Acknowledgements
This work is supported by the [NFDI-Matwerk](https://nfdi-matwerk.de/) consortia.

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under the National Research Data Infrastructure – NFDI 38/1 – project number 460247524

