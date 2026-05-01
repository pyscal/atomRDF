# Changelog

All notable changes to **atomRDF** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — Unreleased

First stable release. The public API documented in [docs/api.rst](docs/api.rst) is
now considered stable; breaking changes will require a major version bump.

If you are upgrading from `0.12.x` (or older `pyscal-rdf`), please read the
[Migration guide](docs/migration.md) — there are several **breaking changes**.

### Added
- `atomrdf.__version__` exposed at the top level (single source of truth in
  `atomrdf/_version.py`).
- Public `WorkflowParser` (`atomrdf.WorkflowParser`) for ingesting external
  workflow descriptions.
- Optional dependency groups so you only install what you use:
  - `pip install "atomrdf[oxigraph]"` — Oxigraph triple store backend
  - `pip install "atomrdf[sqlalchemy]"` — SQLAlchemy-backed store
  - `pip install "atomrdf[materials_project]"` — Materials Project lookups (`mp-api`)
  - `pip install "atomrdf[grainboundary]"` — `aimsgb` + `pymatgen` for grain boundaries
  - `pip install "atomrdf[dislocation]"` — `atomman` for dislocation builders
- CI now tests Python 3.9 — 3.12 (was 3.11 only).

### Changed
- **Breaking — Activity field renaming.** On `Activity` and all subclasses
  (notably `Simulation`):
  - `initial_sample` → `input_sample`
  - `final_sample`   → `output_sample`
- **Breaking — Ontology IRI fix.** The class emitted for crystalline materials
  is now `cdco:CrystallineMaterial` (it is defined in CDCO, not CMSO). Old
  graphs that use `cmso:CrystallineMaterial` will not be queryable with the new
  type. See [Migration guide](docs/migration.md).
- **Breaking — Ontology IRI fix.** `LennardJonesPotential` and `TensileTest`
  now emit the IRIs that actually exist in ASMO
  (`asmo:Lennard-JonesPotential`, `asmo:Tensile_Test`). The Python class names
  are unchanged.
- `Activity` now accepts a `graph=` reference; string sample IDs are
  auto-resolved to `AtomicScaleSample` instances on assignment.
- Library is now silent by default — internal `print()` calls in the GEXF
  exporter were replaced with the standard `logging` module.
- Bare `except:` clauses replaced with typed exceptions plus debug logging in
  `graph.py`, `properties.py`, `build/defect.py`,
  `datamodels/structure.py`, `datamodels/structure_io.py`,
  `datamodels/activity.py`.
- Dependency lists in `setup.py`, `requirements.txt` and `environment.yml`
  reconciled. `networkx` was dropped (unused). Lazy-imported packages
  (`mp-api`, `aimsgb`, `pymatgen`, `atomman`) moved to extras.
- Project description and README rewritten for v1.0.

### Removed
- Stale top-level `index.ipynb` (used a pre-`KnowledgeGraph` API). The book's
  landing page is now `index.md`.
- Dead/commented `get_basis_occupancy` in `atomrdf/properties.py`.
- Unused `import pickle` in `atomrdf/graph.py`.

### Fixed
- `docs/api.rst` no longer references modules that no longer exist
  (`atomrdf.structure`, `atomrdf.network.network`, `atomrdf.workflow.workflow`).
- `docs/license.md` now refers to *atomRDF* (was still saying *pyscal-rdf*).

---

## [0.12.3] and earlier

See the [GitHub release page](https://github.com/pyscal/atomRDF/releases) for
the change history of pre-1.0 versions.
