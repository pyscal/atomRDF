from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='atomrdf',
    version='0.9.1',
    author='Abril Azocar Guzman, Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Ontology based structural manipulation and quering',
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=find_packages(include=['atomrdf', 'atomrdf.*']),
    zip_safe=False,
    download_url = 'https://github.com/pyscal/atomrdf',
    url = 'https://pyscal.org',
    install_requires=['numpy', 'ase', 'rdflib', 
    'pyyaml', 'graphviz', 'networkx', 
    'pyscal3', 'spglib', 'pandas', 'owlready2',
    'atomman', 'mp-api'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
)
