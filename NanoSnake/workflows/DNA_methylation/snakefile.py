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

#~~~~~~~~~~~~~~~~~~~~~~~Load sample sheets~~~~~~~~~~~~~~~~~~~~~~~#
sample_df = pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list = list(sample_df.index)

#~~~~~~~~~~~~~~~~~~~~~~~Create shortcuts~~~~~~~~~~~~~~~~~~~~~~~#
wrappers_dir = config["wrappers_dir"]
merge_fastq_dir = config["merge_fastq"]["outdir"]
fastqc_dir = config["fastqc"]["outdir"]
minimap2_index_dir = config["minimap2_index"]["outdir"]
minimap2_align_dir = config["minimap2_align"]["outdir"]
bamqc_dir = config["bamqc"]["outdir"]
samtools_filter_dir = config["samtools_filter"]["outdir"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def get_fastq (wildcards):
    return glob (sample_df.loc[wildcards.sample, "fastq"])
def get_fast5_dir (wildcards):
    return glob (sample_df.loc[wildcards.sample, "fast5_dir"])
def get_seq_summary (wildcards):
    return glob (sample_df.loc[wildcards.sample, "seq_summary"])

rule all:
    input:
        expand(path.join("results", merge_fastq_dir,"{sample}.fastq"), sample=sample_list),
        expand(path.join("results", fastqc_dir,"{sample}_fastqc.html"), sample=sample_list),
        expand(path.join("results", fastqc_dir,"{sample}_fastqc.zip"), sample=sample_list),
        expand(path.join("results", minimap2_index_dir,"ref.mmi")),
        expand(path.join("results", minimap2_align_dir,"{sample}.bam"), sample=sample_list),
        expand(path.join("results", bamqc_dir,"{sample}","qualimapReport.html"), sample=sample_list),
        expand(path.join("results", bamqc_dir,"{sample}_samtools_stats.txt"), sample=sample_list),
        expand(path.join("results", bamqc_dir,"{sample}_samtools_flagstat.txt"), sample=sample_list),
        expand(path.join("results", bamqc_dir,"{sample}_samtools_idxstats.txt"), sample=sample_list),
        expand(path.join("results", samtools_filter_dir,"{sample}.bam"), sample=sample_list),

rule merge_fastq:
    input:
        get_fastq
    output:
        path.join("results", merge_fastq_dir,"{sample}.fastq")
    log:
        path.join("logs", merge_fastq_dir,"{sample}.log")
    wrapper:
        "file:"+ path.join(wrappers_dir, "concat_files")

rule fastqc:
    input:
        rules.merge_fastq.output
    output:
        html=path.join("results", fastqc_dir,"{sample}_fastqc.html"),
        zip=path.join("results", fastqc_dir,"{sample}_fastqc.zip")
    log:
        path.join("logs", fastqc_dir,"{sample}_fastqc.log")
    params:
        opt=config["fastqc"]["opt"]
    threads:
        config["fastqc"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "fastqc")

rule minimap2_index:
    input:
        config["reference"]
    output:
        path.join("results", minimap2_index_dir,"ref.mmi")
    log:
        path.join("logs", minimap2_index_dir,"ref.log")
    params:
        opt=config["minimap2_index"]["opt"]
    threads:
        config["minimap2_index"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "minimap2_index")

rule minimap2_align:
    input:
        index=rules.minimap2_index.output,
        fastq=rules.merge_fastq.output
    output:
        path.join("results", minimap2_align_dir,"{sample}.bam")
    log:
        path.join("logs", minimap2_align_dir,"{sample}.log")
    params:
        opt=config["minimap2_align"]["opt"],
    threads:
        config["minimap2_align"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "minimap2_align")

rule bamqc:
    input:
        rules.minimap2_align.output,
    output:
        qualimap=path.join("results", bamqc_dir,"{sample}","qualimapReport.html"),
        stats=path.join("results", bamqc_dir,"{sample}_samtools_stats.txt"),
        flagstat=path.join("results", bamqc_dir,"{sample}_samtools_flagstat.txt"),
        idxstats=path.join("results", bamqc_dir,"{sample}_samtools_idxstats.txt"),
    log:
        path.join("logs", bamqc_dir,"{sample}.log")
    wrapper:
        "file:"+ path.join(wrappers_dir, "bamqc")

rule samtools_filter:
    input:
        rules.minimap2_align.output
    output:
        path.join("results", samtools_filter_dir,"{sample}.bam")
    log:
        path.join("logs", samtools_filter_dir,"{sample}.log")
    params:
        opt=config["samtools_filter"]["opt"],
    threads:
        config["samtools_filter"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "samtools_filter")

rule nanopolish_index:
    input:
        fastq = rules.merge_fastq.output
        fast5 = get_fast5_dir
    output:
        fastq.index
    log:
        path.join("logs", samtools_filter_dir,"{sample}.log")
    params:
        opt=config["samtools_filter"]["opt"],
    threads:
        config["samtools_filter"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "samtools_filter")

rule samtools_filter:
    input:
        rules.minimap2_align.output
    output:
        path.join("results", samtools_filter_dir,"{sample}.bam")
    log:
        path.join("logs", samtools_filter_dir,"{sample}.log")
    params:
        opt=config["samtools_filter"]["opt"],
    threads:
        config["samtools_filter"]["threads"]
    wrapper:
        "file:"+ path.join(wrappers_dir, "samtools_filter")
