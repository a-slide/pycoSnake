# Snakemake workflow: NanoSnake

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Licence](https://anaconda.org/aleg/nanosnake/badges/license.svg)](https://anaconda.org/aleg/nanosnake)
[![Build Status](https://travis-ci.com/a-slide/NanoSnake.svg?branch=master)](https://travis-ci.com/a-slide/NanoSnake)
[![DOI](https://zenodo.org/badge/173960745.svg)](https://zenodo.org/badge/latestdoi/173960745)

[![Anaconda version](https://anaconda.org/aleg/nanosnake/badges/version.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda last release](https://anaconda.org/aleg/nanosnake/badges/latest_release_relative_date.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda platforms](https://anaconda.org/aleg/nanosnake/badges/platforms.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda Downloads](https://anaconda.org/aleg/nanosnake/badges/downloads.svg)](https://anaconda.org/aleg/nanosnake)

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


```
conda create -n nanosnake -c aleg -c anaconda -c bioconda -c conda-forge python=3.6 nanosnake

conda activate nanosnake
```

To update the package, from the package directory

```
conda activate nanosnake

conda update nanosnake -c aleg
```

## Usage

At the moment there is only 1 workflow available in NanoSnake:

* DNA : Analyse Basecalled Nanopore sequencing data.
    * Fastq merging and filtering
    * Alignment with minimap2
    * Alignment cleaning
    * Generate coverage plots IGV and bedgraph (optional)
    * Run Nanopore QC with pycoQC (optional)
    * Run methylation analysing with Nanoplolish + pycoMeth (optional)
    * Run SV analysis with Sniffles + pycoSV (optional)

## Configure workflow

Generate the sample sheet and config template files required for the workflow you want to run.

```
conda activate NanoSnake

NanoSnake {WORKFLOW NAME} --generate_template all --overwrite_template

or

NanoSnake {WORKFLOW NAME} -g all -o
```

The `samples.tsv` file needs to be filled with the required informations detailed in the file header and passed to NanoSnake (`--sample_sheet`).

The `config.yaml` can be modified and passed to NanoSnake (`--config`). It is generally recommended to stick to the default parameters.

The `cluster_config.yaml` can be modified and passed to NanoSnake (`--cluster_config`). Use the file instead of config.yaml if you are executing the pipeline in a cluster environment. By default the file is for an LSF cluster, but it can be modified for other HPC platfoms.  

### Step 4: Execute workflow

Call NanoSnake and choose your workflow

```
conda activate NanoSnake

NanoSnake {WORKFLOW NAME} {OPTIONS}
```

#### Example for the DNA workflow

**Usage on a local machine**

```
conda activate NanoSnake

NanoSnake DNA_methylation -r ref.fa  -s sample_sheet.tsv -c config.yaml --cores 10
```

**Usage in an LSF cluster environment**

Use the cluster_config file instead of the normal config

```
conda activate NanoSnake

NanoSnake DNA_methylation  \
    -r ref.fa \
    -s sample_sheet.tsv \
    -c cluster_config.yaml \
    --cores 100 \
    --nodes 5 \
    --cluster "bsub -q {cluster.queue} -n {cluster.threads} -M {cluster.mem} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}"
```
