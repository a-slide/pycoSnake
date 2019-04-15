# Snakemake workflow: NanoSnake

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/a-slide/NanoSnake.svg?branch=master)](https://travis-ci.com/a-slide/NanoSnake#)

---

**UNSTABLE PACKAGE UNDER DEVELOPMENT**

---

## Authors

* Adrien Leger (@a-slide)

## Usage

### Step 1: Install conda package manager

Conda is the only dependency that you need to install manually.

All the other packages and external program needed for this pipeline will then be automatically handled by conda itself.

Install conda following the official documentation for you system

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

### Step 2: Install workflow

#### Develop mode

If you want to install the package in develop mode (installed and editable at the same time), clone the repository locally:

```
git clone git@github.com:a-slide/NanoSnake.git

cd NanoSnake
```

Create a virtual environment to deploy the package

```
conda env create -f environment.yml
```

Activate the conda environment and install NanoSnake with pip

```
conda activate NanoSnake
pip install -e ./
```

#### Update workflow

From the package directory

```
git pull
conda env update -f conda_env.yml
conda activate NanoSnake
pip install -e ./ -U
```

### Step 3: Configure workflow

Generate the sample sheet and config template files required for the workflow you want to run.

```
conda activate NanoSnake

NanoSnake generate_template {PIPELINE NAME}
```

The `samples.tsv` file needs to be filled with the required informations detailed in the file header.

The `config.yaml` can be modified and passed to NanoSnake upon execution, but it is recommended to stick to the default parameters.


### Step 4: Execute workflow

Call NanoSnake and choose your workflow

```
conda activate NanoSnake

NanoSnake {PIPELINE NAME}
```

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
