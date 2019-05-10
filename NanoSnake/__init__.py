# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '0.0.0.a2'
__description__="""NanoSnake contains a collection of snakemake workflows for analysing nanopore sequencing data"""

# Collect info in a dictionnary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __description__,
    "url": "https://github.com/a-slide/NanoSnake",
    "author": 'Adrien Leger',
    "author_email": 'aleg@ebi.ac.uk',
    "license": 'MIT',
    "python_requires":'>=3.6',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'],
    "install_requires": ['pandas>=0.24.1', 'snakemake>=5.4.2', 'tqdm>=4.23'],
    "packages": [__name__], ##############################
    "package_dir": {__name__: __name__}, ##################################
    "package_data": {__name__: ['envs/*', 'rules/*', 'scripts/*', 'snakefiles/*', 'templates/*']}, ##################################
    "entry_points": {'console_scripts': ['NanoSnake=NanoSnake.__main__:main']
    }
}
