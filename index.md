# atomRDF

```{note}
`atomRDF` was previously called `pyscal-rdf`.
```

`atomRDF` is a Python tool for **ontology-based creation, manipulation and querying**
of atomic-scale structures and the simulation workflows that produce them.

The data model is implemented as [Pydantic](https://docs.pydantic.dev/) classes generated from the
[Conceptual Dictionary for Computational Materials Science](https://github.com/OCDO/conceptual_dictionary)
ontologies maintained by [OCDO](https://github.com/OCDO):

| Ontology | What it covers |
| -------- | -------------- |
| **CMSO**  | Computational Material Sample Ontology |
| **CDCO**  | Crystallographic Defect Core Ontology |
| **PODO**  | Point-defect Ontology |
| **PLDO**  | Planar-defect Ontology |
| **LDO**   | Line-defect Ontology |
| **ASMO**  | Atomistic Simulation Methods Ontology |

Every object created by `atomRDF` is round-trippable between Python, JSON/YAML
and RDF, with ontology-conformant semantics out of the box.

## Where to next

- [](docs/gettingstarted.md) — install and run your first example.
- [](docs/examples.md) — worked notebooks (bulk structures, grain boundaries,
  defects, SPARQL queries, working with data, …).
- [](docs/migration.md) — upgrading from `0.12.x`.
- [](docs/api.rst) — full API reference.
- [](docs/extending.md) — contribute / extend.

## Citing atomRDF

If you use atomRDF in academic work, please cite:

> Guzmán, A. A., Menon, S., Hickel, T., & Sandfeld, S. (2026).
> *Ontology-based knowledge graph infrastructure for interoperable atomistic
> simulation data.* arXiv:2604.06230.
> <https://arxiv.org/abs/2604.06230>
