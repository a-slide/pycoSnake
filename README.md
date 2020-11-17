# pycoSnake v0.2.0.dev17

![](pictures/pycoSnake_logo.png)

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Licence](https://anaconda.org/aleg/pycoSnake/badges/license.svg)](https://anaconda.org/aleg/pycoSnake)
[![Build Status](https://travis-ci.com/a-slide/pycoSnake.svg?branch=master)](https://travis-ci.com/a-slide/pycoSnake)
[![DOI](https://zenodo.org/badge/173960745.svg)](https://zenodo.org/badge/latestdoi/173960745)

[![Anaconda version](https://anaconda.org/aleg/pycoSnake/badges/version.svg)](https://anaconda.org/aleg/pycoSnake)
[![Anaconda last release](https://anaconda.org/aleg/pycoSnake/badges/latest_release_relative_date.svg)](https://anaconda.org/aleg/pycoSnake)
[![Anaconda platforms](https://anaconda.org/aleg/pycoSnake/badges/platforms.svg)](https://anaconda.org/aleg/pycoSnake)
[![Anaconda Downloads](https://anaconda.org/aleg/pycoSnake/badges/downloads.svg)](https://anaconda.org/aleg/pycoSnake)
---

**pycoSnake is a neatly wrapped collection of snakemake workflows for analysing Illumina and nanopore sequencing datasets. It is easy to install with conda and simple to run on a local computer or in a cluster environment**

---

## Installation

### Install conda package manager

Conda is the only dependency that you need to install the package.

All the other packages and external program needed for this pipeline will then be automatically handled by conda itself.

Install conda following the official documentation for you system

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html


### Install the package in a conda environment

#### Recommended installation from Anaconda Cloud

You might not need all the extra channels depending on you conda configuration, but you need to add my channel: `aleg`

```
conda create -y -n pycoSnake -c aleg -c anaconda -c bioconda -c conda-forge python=3.6 pycosnake

conda activate pycoSnake
```

To update the package

```
conda activate pycoSnake

conda update pycoSnake -c aleg
```

#### Local installation in develop mode

Clone pycoSnake repository to your local machine

```
git clone git@github.com:a-slide/pycoSnake.git

cd pycoSnake
```

Create a new conda environment

```
conda create -y -n pycoSnake python=3.6

conda activate pycoSnake
```

Install pycoSnake with pip in develop MODE

```
pip install -e ./
```

## pycoSnake workflows

At the moment there is only 2 workflows available in pycoSnake:

* DNA_ONT : Analyse Basecalled Nanopore sequencing data.
    * Download and cleanup genome
    * Fastq merging and filtering
    * Alignment with minimap2
    * Alignment cleaning
    * Generate coverage plots IGV and bedgraph (optional)
    * Run Nanopore QC with pycoQC (optional)
    * Run DNA methylation analysing with Nanopolish + pycoMeth (optional)
    * Run SV analysis with Sniffles + filtering + multi-sample merging (optional)

* RNA_illumina : Analyse Illumina long RNA-seq sequencing data.
    * Download and cleanup genome, transcriptome and annotations
    * Fastq filtering / control / pre-alignment quality control
    * Genome alignment with STAR
    * Summarize STAR counts for all samples
    * Alignment cleaning
    * Transcriptome pseudo-alignment and transcripts quantification with Salmon
    * Summarize Salmon counts for all samples
    * Generate coverage plots IGV and bedgraph (optional)
    * Count reads per genes with featurecounts and summarize counts for all samples (optional)
    * Calculation of FPKM with cufflinks and summarize FPKM for all samples (optional)

### Configure a workflow

Generate the sample sheet and config template files required for the workflow you want to run.

```
conda activate pycoSnake

pycoSnake {WORKFLOW NAME} --generate_template config sample_sheet --overwrite_template

# Or for a cluster environment

pycoSnake {WORKFLOW NAME} --generate_template cluster_config sample_sheet --overwrite_template
```

The `samples.tsv` file needs to be filled with the required informations detailed in the file header and passed to pycoSnake (`--sample_sheet`).

The `config.yaml` can be modified and passed to pycoSnake (`--config`). It is generally recommended to stick to the default parameters.

The `cluster_config.yaml` can be modified and passed to pycoSnake (`--cluster_config`). Use the file instead of config.yaml if you are executing the pipeline in a cluster environment. By default the file is for an LSF cluster, but it can be modified for other HPC platforms.

### Execute a workflow

Call pycoSnake and choose your workflow

```
conda activate pycoSnake

pycoSnake {WORKFLOW NAME} {OPTIONS}
```

#### Example for the DNA workflow

**Usage on a local machine**

```
conda activate pycoSnake

pycoSnake DNA_ONT -r ref.fa  -s sample_sheet.tsv --config config.yaml --cores 10
```

**Usage in an LSF cluster environment**

Use the cluster_config option instead of the config file.
The cluster_config provided with pycoSnake is configured to work on an LSF cluster environment
It contains the "prototype bsub command" to be used by Snakemake  `bsub -q {cluster.queue} -n {cluster.threads} -M {cluster.mem} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}` as well as the maximal number of cores and nodes to use.

```
conda activate pycoSnake

pycoSnake DNA_ONT -r ref.fa -s sample_sheet.tsv --cluster_config cluster_config.yaml
```

## Wrapper library

This repository contains snakemake wrappers for [pycoSnake](https://github.com/a-slide/pycoSnake).

### Using the wrappers library

Wrappers can also be used outside of `pycoSnake` in any Snakemake file by using the following option:

```
--wrapper-prefix https://raw.githubusercontent.com/a-slide/pycoSnake/master/pycoSnake/wrappers/
```

### Testing Wrappers

The package contains test data and integrated tests for all the wrappers.

```
# Test all the wrappers

pycoSnake test_wrappers --cores 8

# Test individual wrappers

pycoSnake test_wrappers --wrappers get_annotation star_count_merge pycometh_comp_report --cores 8

# Keep of test output files generated by the wrappers

pycoSnake test_wrappers --keep_output --cores 8
```

## Command line help message

```
usage: pycoSnake [-h] [--version] {test_wrappers,DNA_ONT,RNA_illumina} ...

pycoSnake is a neatly wrapped collection of snakemake workflows for analysing
nanopore sequencing data

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

subcommands:
  pycoSnake implements the following subcommands

  {test_wrappers,DNA_ONT,RNA_illumina}
```

## Classifiers

* Development Status :: 3 - Alpha
* Intended Audience :: Science/Research
* Topic :: Scientific/Engineering :: Bio-Informatics
* License :: OSI Approved :: MIT License
* Programming Language :: Python :: 3

## citation

Adrien Leger. a-slide/pycoSnake: (2020). doi:10.5281/zenodo.4110611

## licence

MIT

Copyright © 2020 Adrien Leger

## Authors

* Adrien Leger / aleg@ebi.ac.uk / https://adrienleger.com
