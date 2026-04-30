"""atomRDF — ontology-based knowledge graphs for atomistic simulation data.

atomRDF combines `pyscal3 <https://github.com/pyscal/pyscal3>`_,
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ and `RDFLib <https://rdflib.dev>`_ with
the `OCDO Conceptual Dictionary <https://github.com/OCDO/conceptual_dictionary>`_
ontologies (CMSO, CDCO, PODO, PLDO, LDO, ASMO) to make atomic-scale samples and
the simulation workflows that produce them queryable as RDF.

The public entry points re-exported here form the stable v1 API:

* :class:`KnowledgeGraph` -- main container; create samples, add provenance, run
  SPARQL, persist to disk.
* :class:`WorkflowParser` -- ingest external workflow descriptions and add them
  to a :class:`KnowledgeGraph`.

For a quick tour see ``examples/01_getting_started.ipynb`` or the online
documentation at https://atomrdf.pyscal.org.
"""

from atomrdf._version import __version__
from atomrdf.graph import KnowledgeGraph
from atomrdf.io.workflow_parser import WorkflowParser

__all__ = [
    "__version__",
    "KnowledgeGraph",
    "WorkflowParser",
]
