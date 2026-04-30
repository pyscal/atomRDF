# Migration guide

This page explains how to upgrade existing code to **atomRDF v1.0**. Most
users will be coming from `0.12.x` (or older `pyscal-rdf` releases).

If something here is unclear or you find a breaking change that is not
documented, please [open an issue](https://github.com/pyscal/atomRDF/issues).

---

## 1. `Activity` field renaming

The most common breaking change. All `Activity` subclasses (most importantly
`Simulation`) renamed two fields for clarity:

| Before (`<= 0.12`) | After (`>= 1.0`) |
| ------------------ | ----------------- |
| `initial_sample`   | `input_sample`    |
| `final_sample`     | `output_sample`   |

**Before**

```python
sim = Simulation(
    initial_sample=sample_in,
    final_sample=sample_out,
    method=...,
    algorithm=...,
)
```

**After**

```python
sim = Simulation(
    input_sample=sample_in,
    output_sample=sample_out,
    method=...,
    algorithm=...,
)
```

A quick `sed` over your codebase will catch most cases:

```bash
grep -rl "initial_sample\|final_sample" .
```

## 2. Sample IDs are auto-resolved when you pass a graph

You can now hand `Activity` constructors a string sample ID together with a
`graph=` reference and the object will be loaded for you:

```python
sim = Simulation(
    input_sample="sample:1b7006a0-852b-4114-afc8-5b244703500e",
    output_sample="sample:9c8f14a2-963d-4225-b0d2-8e345814611f",
    graph=kg,                  # NEW: enables automatic loading
    method=...,
    algorithm=...,
)
```

Without `graph=`, strings are kept as strings (the old behaviour), so this
change is additive.

## 3. Ontology IRI fixes

Three IRIs emitted to the RDF graph were corrected to match what is actually
published in the upstream ontologies. **The Python class names and the way
you construct objects do not change** — only the type IRI written into the
graph.

| Concept                         | Old IRI (did not exist)               | New IRI (now correct)              |
| ------------------------------- | ------------------------------------- | ---------------------------------- |
| Crystalline material            | `cmso:CrystallineMaterial`            | `cdco:CrystallineMaterial`         |
| Lennard-Jones potential         | `asmo:LennardJonesPotential`          | `asmo:Lennard-JonesPotential`      |
| Tensile-test algorithm          | `asmo:TensileTest`                    | `asmo:Tensile_Test`                |

### Impact on existing graphs

Persisted graphs (`.ttl`, Oxigraph store, SQLAlchemy DB, …) created with
`<= 0.12` use the *old* IRIs. SPARQL queries that filter by the new IRI will
not return those triples and vice-versa. You have two options:

**Option A — re-emit the graph from your source data** (recommended if your
data is small or regeneratable).

**Option B — patch the graph in place with a SPARQL UPDATE:**

```sparql
PREFIX cmso:  <http://purls.helmholtz-metadaten.de/cmso/>
PREFIX cdco:  <http://purls.helmholtz-metadaten.de/cdos/cdco/>
PREFIX asmo:  <http://purls.helmholtz-metadaten.de/asmo/>
PREFIX rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

DELETE { ?s rdf:type cmso:CrystallineMaterial }
INSERT { ?s rdf:type cdco:CrystallineMaterial }
WHERE  { ?s rdf:type cmso:CrystallineMaterial };

DELETE { ?s rdf:type asmo:LennardJonesPotential }
INSERT { ?s rdf:type <http://purls.helmholtz-metadaten.de/asmo/Lennard-JonesPotential> }
WHERE  { ?s rdf:type asmo:LennardJonesPotential };

DELETE { ?s rdf:type asmo:TensileTest }
INSERT { ?s rdf:type asmo:Tensile_Test }
WHERE  { ?s rdf:type asmo:TensileTest };
```

## 4. Optional dependencies are now extras

Heavy or rarely-used dependencies are no longer installed by default. If you
use these features, install the matching extra:

| Feature                                              | Install                                  |
| ---------------------------------------------------- | ---------------------------------------- |
| Materials Project lookups (`mp-api`)                 | `pip install "atomrdf[materials_project]"` |
| Grain boundary builders (`aimsgb`, `pymatgen`)       | `pip install "atomrdf[grainboundary]"`     |
| Dislocation builders (`atomman`)                     | `pip install "atomrdf[dislocation]"`       |
| Oxigraph triple store backend                        | `pip install "atomrdf[oxigraph]"`          |
| SQLAlchemy-backed store                              | `pip install "atomrdf[sqlalchemy]"`        |

Combine with commas: `pip install "atomrdf[grainboundary,oxigraph]"`.

## 5. Library is silent by default

The GEXF exporter no longer prints to stdout. If you want the old summary,
configure logging:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

## 6. Removed / renamed bits

- The legacy top-level `index.ipynb` (which used the pre-`KnowledgeGraph` API)
  has been removed.
- `atomrdf.properties.get_basis_occupancy` (long commented out) has been
  removed for good.
- `atomrdf.graph` no longer imports `pickle` at the top level.

## 7. Python support

Python `>=3.9, <3.13` is now tested in CI. Earlier versions may still work
but are no longer supported.
