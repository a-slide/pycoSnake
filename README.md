# Snakemake workflow: NanoSnake

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/NanoSnake.svg?branch=master)](https://travis-ci.org/snakemake-workflows/NanoSnake)

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

#### Initial installation

If you simply want to use this workflow, clone the repository locally:

`git clone git@github.com:a-slide/NanoSnake.git`

`cd NanoSnake`

Then create a virtual environment to deploy the workflow

`conda env create -f conda_env.yml`


#### Update workflow

From the package directory

`git pull`

`conda env update -f conda_env.yml`

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

First activate the conda virtual environment

`conda activate NanoSnake`

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N


See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.
