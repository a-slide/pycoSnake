![](pictures/NanoSnake_logo.png)

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Licence](https://anaconda.org/aleg/nanosnake/badges/license.svg)](https://anaconda.org/aleg/nanosnake)
[![Build Status](https://travis-ci.com/a-slide/NanoSnake.svg?branch=master)](https://travis-ci.com/a-slide/NanoSnake)
[![DOI](https://zenodo.org/badge/173960745.svg)](https://zenodo.org/badge/latestdoi/173960745)

[![Anaconda version](https://anaconda.org/aleg/nanosnake/badges/version.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda last release](https://anaconda.org/aleg/nanosnake/badges/latest_release_relative_date.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda platforms](https://anaconda.org/aleg/nanosnake/badges/platforms.svg)](https://anaconda.org/aleg/nanosnake)
[![Anaconda Downloads](https://anaconda.org/aleg/nanosnake/badges/downloads.svg)](https://anaconda.org/aleg/nanosnake)

---

**NanoSnake is a collection of [Snakemake](https://github.com/snakemake/snakemake) pipelines neatly wrapped in a convenient python package interface.
It is easy to install with conda and to simple to run on a local computer or in a cluster environment.**

---

## Installation

### Install conda package manager

Conda is the only dependency that you need to install the package.

All the other packages and external program needed for this pipeline will then be automatically handled by conda itself.

Install conda following the official documentation for you system

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

### Install the package in a conda environment

You might not need all the extra channels depending on you conda configuration, but you need to add at least my channel: `aleg`

```
conda create -n nanosnake -c aleg -c anaconda -c bioconda -c conda-forge python=3.6 nanosnake

conda activate nanosnake
```

To update the package, from the package directory

```
conda activate nanosnake

conda update nanosnake -c aleg
```

## Running NanoSnake workflows

At the moment there is only 2 workflows available in NanoSnake:

* DNA_ONT : Analyse Basecalled Nanopore sequencing data.
    * Download and cleanup genome
    * Fastq merging and filtering
    * Alignment with minimap2
    * Alignment cleaning
    * Generate coverage plots IGV and bedgraph (optional)
    * Run Nanopore QC with pycoQC (optional)
    * Run DNA methylation analysing with Nanoplolish + pycoMeth (optional)
    * Run SV analysis with Sniffles + pycoSV (optional)

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
conda activate NanoSnake

NanoSnake {WORKFLOW NAME} --generate_template config sample_sheet --overwrite_template

# Or for a cluster environment

NanoSnake {WORKFLOW NAME} --generate_template cluster_config sample_sheet --overwrite_template
```

The `samples.tsv` file needs to be filled with the required informations detailed in the file header and passed to NanoSnake (`--sample_sheet`).

The `config.yaml` can be modified and passed to NanoSnake (`--config`). It is generally recommended to stick to the default parameters.

The `cluster_config.yaml` can be modified and passed to NanoSnake (`--cluster_config`). Use the file instead of config.yaml if you are executing the pipeline in a cluster environment. By default the file is for an LSF cluster, but it can be modified for other HPC platforms.

### Execute a workflow

Call NanoSnake and choose your workflow

```
conda activate NanoSnake

NanoSnake {WORKFLOW NAME} {OPTIONS}
```

#### Example for the DNA workflow

**Usage on a local machine**

```
conda activate NanoSnake

NanoSnake DNA_ONT -r ref.fa  -s sample_sheet.tsv --config config.yaml --cores 10
```

**Usage in an LSF cluster environment**

Use the cluster_config option instead of the config file.
The cluster_config provided with NanoSnake is configured to work on an LSF cluster environment
It contains the "prototype bsub command" to be used by Snakemake  `bsub -q {cluster.queue} -n {cluster.threads} -M {cluster.mem} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}` as well as the maximal number of cores and nodes to use.

```
conda activate NanoSnake

NanoSnake DNA_ONT -r ref.fa -s sample_sheet.tsv --cluster_config cluster_config.yaml
```

## Wrapper library

This repository contains snakemake wrappers for [NanoSnake](https://github.com/a-slide/NanoSnake).

### Using the wrappers library

Wrappers can also be used outside of `NanoSnake` in any Snakemake file by using the following option:

```
--wrapper-prefix https://raw.githubusercontent.com/a-slide/NanoSnake/master/NanoSnake/wrappers/
```

### Testing Wrappers

The package contains test data and integrated tests all the wrappers.

### Command line help message

```
usage: NanoSnake tests [-h]
                       [--wrappers ...]]
                       [--keep_output] [--clean_output] [--cores CORES]
                       [--workdir WORKDIR] [--verbose | --quiet]

Test all the wrappers

optional arguments:
  -h, --help                    show this help message and exit
  --wrappers ...]               List of wrappers to test (default: all)
  --keep_output, -k             Keep temporary output files generated during tests (default: False)
  --clean_output, -c            clean all temporary output files generated during tests (default: False)
  --cores CORES, -j CORES       the number of provided cores (default: 1)
  --workdir WORKDIR, -d WORKDIR Path to the working dir where to deploy the workflow (default: ./)
  --verbose, -v                 Show additional debug output (default: False)
  --quiet, -q                   Reduce overall output (default: False)
```

### Citing

The repository is archived at Zenodo. If you use `NanoSnake` please cite as follow:

Adrien Leger. (2019, July 30). a-slide/NanoSnake. Zenodo. http://doi.org/10.5281/zenodo.3352283

## Authors and contributors

* Adrien Leger (@a-slide) - aleg {at} ebi.ac.uk
