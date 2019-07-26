#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name= "NanoSnake"
version = "0.0.0.4"
description = "NanoSnake is a neatly wrapped collection of snakemake workflows for analysing nanopore sequencing data"
with open("README.md", "r") as fh:
    long_description = fh.read()

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
    install_requires = ['pandas>=0.24.1', 'snakemake>=5.5.4'],
    packages = [name],
    package_dir = {name: name},
    package_data = {name: ['workflows/*', 'wrappers/*']},
    entry_points = {'console_scripts': ['NanoSnake=NanoSnake.__main__:main']}
)
