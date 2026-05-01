from setuptools import setup, find_packages
import os
import re


with open("README.md") as readme_file:
    readme = readme_file.read()


def _read_version():
    """Read ``__version__`` from atomrdf/_version.py without importing the package."""
    here = os.path.abspath(os.path.dirname(__file__))
    version_file = os.path.join(here, "atomrdf", "_version.py")
    with open(version_file, "r") as fh:
        match = re.search(r'^__version__\s*=\s*"([^"]+)"', fh.read(), re.M)
    if not match:
        raise RuntimeError("Unable to find __version__ in atomrdf/_version.py")
    return match.group(1)


setup(
    name="atomrdf",
    version=_read_version(),
    author="Abril Azocar Guzman, Sarath Menon",
    author_email="sarath.menon@pyscal.org",
    description="Ontology-based creation, manipulation and querying of atomic-scale structures and simulation workflows, with Pydantic data models generated from the OCDO Conceptual Dictionary (CMSO, CDCO, PODO, PLDO, LDO, ASMO).",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=find_packages(include=["atomrdf", "atomrdf.*"]),
    zip_safe=False,
    download_url="https://github.com/pyscal/atomrdf",
    url="https://pyscal.org",
    install_requires=[
        # Core scientific stack
        "numpy",
        "pandas",
        # Atomic structures
        "ase",
        "pyscal3",
        "spglib",
        "mendeleev",
        # RDF / ontology tooling
        "rdflib",
        "tools4rdf",
        # Data models / IO
        "pydantic",
        "pyyaml",
        # Visualisation
        "graphviz",
    ],
    extras_require={
        # Triple stores
        "oxigraph": ["oxrdflib"],
        "sqlalchemy": ["rdflib-sqlalchemy"],
        # Optional builders / data sources (lazy-imported)
        "materials_project": ["mp-api"],
        "grainboundary": ["aimsgb", "pymatgen"],
        "dislocation": ["atomman"],
    },
    python_requires=">=3.9",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    include_package_data=True,
)
