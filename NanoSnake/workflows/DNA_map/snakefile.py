# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Std lib
from glob import glob
from os import path
import yaml

# Third party lib
import pandas as pd
from snakemake.utils import min_version

# set minimum snakemake version
min_version("5.4.2")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Load sample sheets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sample_df = pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list = list(sample_df.index)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def get_fastq (wildcards):
    return sample_df.loc[wildcards.sample, "fastq"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Top Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input:
        expand(path.join("results","merge_filter_fastq","{sample}.fastq"), sample=sample_list),
        expand(path.join("results","fastqc","{sample}_fastqc.html"), sample=sample_list),
        expand(path.join("results","fastqc","{sample}_fastqc.zip"), sample=sample_list),
        expand(path.join("results","minimap2_index","ref.mmi")),
        expand(path.join("results","minimap2_align", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","bamqc","{sample}","qualimapReport.html"), sample=sample_list),
        expand(path.join("results","bamqc","{sample}_samtools_stats.txt"), sample=sample_list),
        expand(path.join("results","bamqc","{sample}_samtools_flagstat.txt"), sample=sample_list),
        expand(path.join("results","bamqc","{sample}_samtools_idxstats.txt"), sample=sample_list),
        expand(path.join("results","samtools_filter", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","genomecov","{sample}.bedgraph"), sample=sample_list),

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule merge_filter_fastq:
    input: get_fastq
    output: path.join("results","merge_filter_fastq","{sample}.fastq")
    log: path.join("logs","merge_filter_fastq","{sample}.log")
    threads: config["merge_filter_fastq"].get("threads", 1)
    params: opt=config["merge_filter_fastq"].get("opt", "")
    resources: mem_mb=config["merge_filter_fastq"].get("mem", 1000)
    wrapper: "merge_filter_fastq"

rule fastqc:
    input: rules.merge_filter_fastq.output
    output:
        html=path.join("results","fastqc","{sample}_fastqc.html"),
        zip=path.join("results","fastqc","{sample}_fastqc.zip")
    log: path.join("logs", "fastqc","{sample}_fastqc.log")
    threads: config["fastqc"].get("threads", 2)
    params: opt=config["fastqc"].get("opt", "")
    resources: mem_mb=config["fastqc"].get("mem", 1000)
    wrapper: "fastqc"

rule minimap2_index:
    input: config["reference"]
    output: path.join("results","minimap2_index","ref.mmi")
    log: path.join("logs","minimap2_index","ref.log")
    threads: config["minimap2_index"].get("threads", 1)
    params: opt=config["minimap2_index"].get("opt", "")
    resources: mem_mb=config["minimap2_index"].get("mem", 1000)
    wrapper: "minimap2_index"

rule minimap2_align:
    input:
        index=rules.minimap2_index.output,
        fastq=rules.merge_filter_fastq.output
    output: path.join("results","minimap2_align","{sample}.bam")
    log: path.join("logs","minimap2_align","{sample}.log")
    threads: config["minimap2_align"].get("threads", 4)
    params: opt=config["minimap2_align"].get("opt", "")
    resources: mem_mb=config["minimap2_align"].get("mem", 1000)
    wrapper: "minimap2_align"

rule bamqc:
    input: rules.minimap2_align.output,
    output:
        qualimap=path.join("results","bamqc","{sample}","qualimapReport.html"),
        stats=path.join("results","bamqc","{sample}_samtools_stats.txt"),
        flagstat=path.join("results","bamqc","{sample}_samtools_flagstat.txt"),
        idxstats=path.join("results","bamqc","{sample}_samtools_idxstats.txt"),
    log: path.join("logs","bamqc","{sample}.log")
    threads: config["bamqc"].get("threads", 1)
    params: opt=config["bamqc"].get("opt", "")
    resources: mem_mb=config["bamqc"].get("mem", 1000)
    wrapper: "bamqc"

rule samtools_filter:
    input: rules.minimap2_align.output
    output: path.join("results","samtools_filter","{sample}.bam")
    log: path.join("logs","samtools_filter","{sample}.log")
    threads: config["samtools_filter"].get("threads", 2)
    params: opt=config["samtools_filter"].get("opt", "")
    resources: mem_mb=config["samtools_filter"].get("mem", 1000)
    wrapper: "samtools_filter"

rule genomecov:
    input: rules.samtools_filter.output
    output: path.join("results","genomecov","{sample}.bedgraph")
    log: path.join("logs","genomecov","{sample}.log")
    threads: config["genomecov"].get("threads", 1)
    params: opt=config["genomecov"].get("opt", "")
    resources: mem_mb=config["genomecov"].get("mem", 1000)
    wrapper: "genomecov"
