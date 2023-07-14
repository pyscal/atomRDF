from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='pyscal_rdf',
    version='0.0.14',
    author='Abril Azocar Guzman, Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Ontology based structural manipulation and quering',
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=find_packages(include=['pyscal_rdf', 'pyscal_rdf.*']),
    zip_safe=False,
    download_url = 'https://github.com/pyscal/pyscal_rdf',
    url = 'https://pyscal.org',
    install_requires=['numpy', 'ase', 'rdflib', 
    'pyyaml', 'graphviz', 'networkx', 
    'ipycytoscape', 'pyscal3', 'spglib',],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
)
