#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from glob import glob

# Define package info
name= "NanoSnake"
version = '0.0.3.10'
description = "NanoSnake is a neatly wrapped collection of snakemake workflows for analysing nanopore sequencing data"
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect all package data to add
package_data = []
package_data.extend([fn.partition("/")[-1] for fn in glob("NanoSnake/test_data/**", recursive=True)])
package_data.extend([fn.partition("/")[-1] for fn in glob("NanoSnake/workflows/**", recursive=True)])
package_data.extend([fn.partition("/")[-1] for fn in glob("NanoSnake/wrappers/**", recursive=True)])

setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/a-slide/NanoSnake",
    author = 'Adrien Leger',
    author_email = 'aleg@ebi.ac.uk',
    license = 'MIT',
    python_requires ='>=3.6',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'pandas==0.24.1',
        'snakemake==5.5.4',
        "loguru==0.3.2",
        "ftputil==3.4"],
    packages = [name],
    package_dir = {name: name},
    package_data = {name: package_data},
    entry_points = {'console_scripts': ['NanoSnake=NanoSnake.__main__:main']}
)
