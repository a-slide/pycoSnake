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
ENVS_DIR = config["envs_dir"]
SCRIPTS_DIR = config["scripts_dir"]
RULES_DIR = config["rules_dir"]

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
        expand(path.join("results", config["fastQC"]["outdir"],"{sample}.html"), sample=SAMPLES),
        expand(path.join("results", config["minimap2_index"]["outdir"],"ref.mmi")),
        expand(path.join("results", config["minimap2_align"]["outdir"],"{sample}.sam"), sample=SAMPLES),
        expand(path.join("results", config["samtools_stats"]["outdir"],"{sample}.txt"), sample=SAMPLES),

rule merge_fastq:
    input:
        get_fastq
    output:
        path.join("results", config["merge_fastq"]["outdir"],"{sample}.fastq.gz")
    log:
        path.join("logs", config["merge_fastq"]["outdir"],"{sample}.log")
    shell:
        "cat {input} | gzip -c > {output} 2> {log}"

rule fastQC:
    input:
        rules.merge_fastq.output
    output:
        html=path.join("results", config["fastQC"]["outdir"],"{sample}.html"),
        zip=path.join("results", config["fastQC"]["outdir"],"{sample}.zip")
    params:
        config["fastQC"]["params"]
    log:
        path.join("logs", config["fastQC"]["outdir"],"{sample}.log")
    wrapper:
        "0.31.1/bio/fastqc"

rule minimap2_index:
    input:
        target=config["reference"]
    output:
        path.join("results", config["minimap2_index"]["outdir"],"ref.mmi")
    log:
        path.join("logs", config["minimap2_index"]["outdir"],"index.log")
    params:
        extra=config["minimap2_index"]["params"]
    #threads: 3
    wrapper:
        "0.31.1/bio/minimap2/index"

rule minimap2_align:
    input:
        target=rules.minimap2_index.output,
        query= rules.merge_fastq.output
    output:
        path.join("results", config["minimap2_align"]["outdir"],"{sample}.sam")
    log:
        path.join("logs", config["minimap2_align"]["outdir"],"{sample}.log")
    params:
        extra=config["minimap2_align"]["params"]
    wrapper:
        "0.31.1/bio/minimap2/aligner"

rule samtools_stats:
    input:
        rules.minimap2_align.output,
    output:
        path.join("results", config["samtools_stats"]["outdir"],"{sample}.txt")
    log:
        path.join("logs", config["samtools_stats"]["outdir"],"{sample}.log")
    params:
        extra=config["samtools_stats"]["params"]
    wrapper:
        "0.31.1/bio/samtools/stats"
