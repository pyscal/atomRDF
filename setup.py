from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='pyscal_rdf',
    version='0.0.0-prod0',
    author='Abril Azocar Guzman, Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Ontology based structural manipulation and quering',
    long_description=readme,
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('pyscal_rdf'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'pyscal_rdf'},
    zip_safe=False,
    download_url = 'https://github.com/pyscal/pyscal_rdf',
    url = 'https://pyscal.org',
    install_requires=['numpy', 'ase', 'rdflib'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
    #package_data={'': ['*.yaml']},
)
