#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from glob import glob

# Long description from README file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect all package data to add
package_data = []
package_data.extend([fn.partition("/")[-1] for fn in glob("pycoSnake/test_data/**", recursive=True)])
package_data.extend([fn.partition("/")[-1] for fn in glob("pycoSnake/workflows/**", recursive=True)])
package_data.extend([fn.partition("/")[-1] for fn in glob("pycoSnake/wrappers/**", recursive=True)])

# Collect info in a dictionary for setup.py
setup(
    name="pycoSnake",
    description="pycoSnake is a neatly wrapped collection of snakemake workflows for analysing Illumina and nanopore sequencing datasets. It is easy to install with conda and simple to run on a local computer or in a cluster environment",
    version="0.2.5.post1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/a-slide/pycoSnake",
    author="Adrien Leger",
    author_email="aleg@ebi.ac.uk",
    license="MIT",
    python_requires=">=3.6",
    classifiers=["Development Status :: 3 - Alpha", "Intended Audience :: Science/Research", "Topic :: Scientific/Engineering :: Bio-Informatics", "License :: OSI Approved :: MIT License", "Programming Language :: Python :: 3"],
    install_requires=["snakemake", "pandas", "ftputil"],
    packages=["pycoSnake"],
    package_dir = {"pycoSnake": "pycoSnake"},
    package_data = {"pycoSnake": package_data},
    entry_points={"console_scripts": ["pycoSnake=pycoSnake.__main__:main"]},
)
