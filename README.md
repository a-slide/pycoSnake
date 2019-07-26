# Snakemake workflow: NanoSnake

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Anaconda version](https://anaconda.org/aleg/nanosnake/badges/version.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda last release](https://anaconda.org/aleg/nanosnake/badges/latest_release_relative_date.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda platforms](https://anaconda.org/aleg/nanosnake/badges/platforms.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda Downloads](https://anaconda.org/aleg/nanosnake/badges/downloads.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda Licence](https://anaconda.org/aleg/nanosnake/badges/license.svg)](https://anaconda.org/aleg/nanosnake)

---

**UNSTABLE PACKAGE UNDER ACTIVE DEVELOPMENT**

---

## Authors

* Adrien Leger (@a-slide)

## Installation

### Install conda package manager

Conda is the only dependency that you need to install the package.

All the other packages and external program needed for this pipeline will then be automatically handled by conda itself.

Install conda following the official documentation for you system

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

### Install the package in a conda environment

If you want to install the package in develop mode (installed and editable at the same time), clone the repository locally:

```
```


To update the package, from the package directory

```
```

## Usage

At the moment there are 2 workflows available in NanoSnake:
* DNA_map
* DNA_methylation

## Configure workflow

Generate the sample sheet and config template files required for the workflow you want to run.

```
conda activate NanoSnake

NanoSnake {WORKFLOW NAME} --generate_template --overwrite_template
```

The `samples.tsv` file needs to be filled with the required informations detailed in the file header.

The `config.yaml` can be modified and passed to NanoSnake upon execution, but it is recommended to stick to the default parameters.


### Step 4: Execute workflow

Call NanoSnake and choose your workflow

```
conda activate NanoSnake

NanoSnake {PIPELINE NAME} {OPTIONS}
```
