# pyscal_rdf

`pyscal_rdf` is a python tool for ontology-based creation, manipulation, and quering of structures. `pyscal_rdf` uses the [Computational Material Sample Ontology (CMSO)](https://github.com/Materials-Data-Science-and-Informatics/cmso-ontology). 

The package is currently under activate development and could be unstable.

You can try `pyscal_rdf` here:

| Jupyter notebook  | GUI |
|-------------------|-----|
| [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pyscal/pyscal_rdf/HEAD?labpath=example.ipynb)  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pyscal/pyscal_rdf/voila?urlpath=voila%2Frender%2Fexample_gui.ipynb)  |

## Installation

### Supported operating systems

`pyscal_rdf` can be installed on Linux and Mac OS based systems. On Windows systems, it is recommended to use  [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

### Using pip

```
pip install pyscal-rdf
```
### Using conda

```
conda install -c conda-forge pyscal-rdf
```


### Building from the repository

We strongly recommend creating a conda environment for the installation. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, the following steps will help set up an environment to use `pyscal_rdf`. First step is to clone the repository.

```
git clone https://github.com/pyscal/pyscal_rdf.git
```

After cloning, an environment can be created from the included file-

```
cd pyscal_rdf
conda env create -f environment.yml
```

This will install the necessary packages and create an environment called rdf. It can be activated by,

```
conda activate rdf
```

then, install `pyscal_rdf` using,

```
pip install .
```

## Using `pyscal_rdf`

Coming soon..


## Acknowledgements

This work is supported by the [NFDI-Matwerk](https://nfdi-matwerk.de/) consortia.
