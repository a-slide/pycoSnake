# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Std lib
from glob import glob
from os import path

# Third party lib
import pandas as pd
from snakemake.utils import min_version

# set minimum snakemake version
min_version("5.4.2")

#~~~~~~~~~~~~~~~~~~~~~~~Load sample sheets and shortcuts~~~~~~~~~~~~~~~~~~~~~~~#

sample_df = pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
SAMPLES = sample_df.index
WRAPPERS_DIR = config["WRAPPERS_DIR"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def get_fastq (wildcards):
    return glob (sample_df.loc[wildcards.sample, "fastq"])
def get_fast5_dir (wildcards):
    return glob (sample_df.loc[wildcards.sample, "fast5_dir"])
def get_seq_summary (wildcards):
    return glob (sample_df.loc[wildcards.sample, "seq_summary"])

rule all:
    input:
        expand(path.join("results", config["merge_fastq"]["outdir"],"{sample}.fastq.gz"), sample=SAMPLES),
        expand(path.join("results", config["fastqc"]["outdir"],"{sample}_fastqc.html"), sample=SAMPLES),
        expand(path.join("results", config["fastqc"]["outdir"],"{sample}_fastqc.zip"), sample=SAMPLES),
        expand(path.join("results", config["minimap2_index"]["outdir"],"ref.mmi")),
        expand(path.join("results", config["minimap2_align"]["outdir"],"{sample}.bam"), sample=SAMPLES),
        # expand(path.join("results", config["samtools_stats"]["outdir"],"{sample}.txt"), sample=SAMPLES),

rule concat_files:
    input:
        get_fastq
    output:
        path.join("results", config["merge_fastq"]["outdir"],"{sample}.fastq.gz")
    log:
        path.join("logs", config["merge_fastq"]["outdir"],"{sample}.log")
    wrapper:
        "file:"+ path.join(WRAPPERS_DIR, "concat_files")

rule fastqc:
    input:
        rules.concat_files.output
    output:
        html=path.join("results", config["fastqc"]["outdir"],"{sample}_fastqc.html"),
        zip=path.join("results", config["fastqc"]["outdir"],"{sample}_fastqc.zip")
    log:
        path.join("results", config["fastqc"]["outdir"],"{sample}_fastqc.log")
    params:
        opt=config["fastqc"]["opt"]
    threads:
        config["fastqc"]["threads"]
    wrapper:
        "file:"+ path.join(WRAPPERS_DIR, "fastqc")

rule minimap2_index:
    input:
        config["reference"]
    output:
        path.join("results", config["minimap2_index"]["outdir"],"ref.mmi")
    log:
        path.join("logs", config["minimap2_index"]["outdir"],"ref.log")
    params:
        opt=config["minimap2_index"]["opt"]
    threads:
        config["minimap2_index"]["threads"]
    wrapper:
        "file:"+ path.join(WRAPPERS_DIR, "minimap2", "index")

rule minimap2_align:
    input:
        index=rules.minimap2_index.output,
        fastq=rules.concat_files.output
    output:
        path.join("results", config["minimap2_align"]["outdir"],"{sample}.bam")
    log:
        path.join("logs", config["minimap2_align"]["outdir"],"{sample}.log")
    params:
        opt=config["minimap2_align"]["opt"],
    threads:
        config["minimap2_align"]["threads"]
    wrapper:
        "file:"+ path.join(WRAPPERS_DIR, "minimap2", "align")

rule samtools_qc:
    input:
        rules.minimap2_align.output.bam,
    output:
        stats=path.join("results", config["samtools_qc"]["outdir"],"{sample}_stats.txt"),
        flagstat=path.join("results", config["samtools_qc"]["outdir"],"{sample}_flagstat.txt"),
        idxstats=path.join("results", config["samtools_qc"]["outdir"],"{sample}_idxstats.txt"),
    log:
        path.join("logs", config["samtools_qc"]["outdir"],"{sample}.log")
    params:
        opt=config["samtools_qc"]["opt"]
    wrapper:
        "file:"+ path.join(WRAPPERS_DIR, "samtools_qc")
