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

def get_fast5 (wildcards):
    return sample_df.loc[wildcards.sample, "fast5"]

def get_seq_summary (wildcards):
    return sample_df.loc[wildcards.sample, "seq_summary"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Top Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input:
        expand(path.join("results","pbt_fastq_filter","{sample}.fastq"), sample=sample_list),
        expand(path.join("results","minimap2_index","ref.mmi")),
        expand(path.join("results","minimap2_align", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","pbt_alignment_filter", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","pbt_fastq_filter","{sample}.fastq.index"), sample=sample_list) if "nanopolish_index" in config else [],
        expand(path.join("results", "nanopolish_call_methylation","{sample}.tsv"), sample=sample_list) if "nanopolish_call_methylation" in config else [],
        expand(path.join("results", "pycometh_aggregate","{sample}.bed"), sample=sample_list) if "pycometh_aggregate" in config else [],
        expand(path.join("results", "pycometh_aggregate","{sample}.tsv"), sample=sample_list) if "pycometh_aggregate" in config else [],
        expand(path.join("results", "sniffles","{sample}.vcf"), sample=sample_list) if "sniffles" in config else [],
        expand(path.join("results","pycoqc","{sample}_pycoqc.html"), sample=sample_list) if "pycoqc" in config else [],
        expand(path.join("results","pycoqc","{sample}_pycoqc.json"), sample=sample_list) if "pycoqc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_stats.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_flagstat.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_idxstats.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","genomecov","{sample}.bedgraph"), sample=sample_list) if "genomecov" in config else [],
        expand(path.join("results","igvtools_count","{sample}.tdf"), sample=sample_list) if "igvtools_count" in config else [],

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CORE RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule pbt_fastq_filter:
    input: get_fastq
    output: path.join("results","pbt_fastq_filter","{sample}.fastq")
    log: path.join("logs","pbt_fastq_filter","{sample}.log")
    threads: config["pbt_fastq_filter"].get("threads", 1)
    params: opt=config["pbt_fastq_filter"].get("opt", "")
    resources: mem_mb=config["pbt_fastq_filter"].get("mem", 1000)
    wrapper: "pbt_fastq_filter"

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
        fastq=rules.pbt_fastq_filter.output
    output: path.join("results","minimap2_align","{sample}.bam")
    log: path.join("logs","minimap2_align","{sample}.log")
    threads: config["minimap2_align"].get("threads", 4)
    params: opt=config["minimap2_align"].get("opt", "")
    resources: mem_mb=config["minimap2_align"].get("mem", 1000)
    wrapper: "minimap2_align"

rule pbt_alignment_filter:
    input: rules.minimap2_align.output
    output: path.join("results","pbt_alignment_filter","{sample}.bam")
    log: path.join("logs","pbt_alignment_filter","{sample}.log")
    threads: config["pbt_alignment_filter"].get("threads", 1)
    params: opt=config["pbt_alignment_filter"].get("opt", "")
    resources: mem_mb=config["pbt_alignment_filter"].get("mem", 1000)
    wrapper: "pbt_alignment_filter"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DNA METHYLATION RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "nanopolish_index" in config:
    rule nanopolish_index:
        input:
            fastq = rules.pbt_fastq_filter.output,
            fast5 = get_fast5,
            seq_summary = get_seq_summary,
        output: path.join("results","pbt_fastq_filter","{sample}.fastq.index")
        log: path.join("logs","nanopolish_index","{sample}.log")
        threads: config["nanopolish_index"].get("threads", 2)
        params: opt=config["nanopolish_index"].get("opt", ""),
        resources: mem_mb=config["nanopolish_index"].get("mem", 1000)
        wrapper: "nanopolish_index"

if "nanopolish_call_methylation" in config:
    rule nanopolish_call_methylation:
        input:
            fastq = rules.pbt_fastq_filter.output,
            fastq_index = rules.nanopolish_index.output,
            bam = rules.pbt_alignment_filter.output,
            ref = config["reference"],
        output: path.join("results", "nanopolish_call_methylation","{sample}.tsv")
        log: path.join("logs","nanopolish_call_methylation","{sample}.log")
        threads: config["nanopolish_call_methylation"].get("threads", 2)
        params: opt=config["nanopolish_call_methylation"].get("opt", ""),
        resources: mem_mb=config["nanopolish_call_methylation"].get("mem", 1000)
        wrapper: "nanopolish_call_methylation"

if "pycometh_aggregate" in config:
    rule pycometh_aggregate:
        input:
            call = rules.nanopolish_call_methylation.output,
            ref = config["reference"],
        output:
            bed = path.join("results", "pycometh_aggregate","{sample}.bed"),
            tsv = path.join("results", "pycometh_aggregate","{sample}.tsv"),
        log: path.join("logs","pycometh_aggregate","{sample}.log")
        threads: config["pycometh_aggregate"].get("threads", 1)
        params: opt=config["pycometh_aggregate"].get("opt", ""),
        resources: mem_mb=config["pycometh_aggregate"].get("mem", 1000)
        wrapper: "pycometh_aggregate"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QC RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "sniffles" in config:
    rule sniffles:
        input: rules.minimap2_align.output
        output: path.join("results", "sniffles","{sample}.vcf")
        log: path.join("logs","sniffles","{sample}.log")
        threads: config["sniffles"].get("threads", 1)
        params: opt=config["sniffles"].get("opt", ""),
        resources: mem_mb=config["sniffles"].get("mem", 1000)
        wrapper: "sniffles"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QC RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "pycoqc" in config:
    rule pycoqc:
        input:
            seq_summary = get_seq_summary,
            bam = rules.minimap2_align.output
        output:
            html=path.join("results","pycoqc","{sample}_pycoqc.html"),
            json=path.join("results","pycoqc","{sample}_pycoqc.json")
        log: path.join("logs", "pycoqc","{sample}_pycoqc.log")
        threads: config["pycoqc"].get("threads", 2)
        params: opt=config["pycoqc"].get("opt", "")
        resources: mem_mb=config["pycoqc"].get("mem", 1000)
        wrapper: "pycoqc"

if "samtools_qc" in config:
    rule samtools_qc:
        input: rules.minimap2_align.output,
        output:
            stats=path.join("results","samtools_qc","{sample}_samtools_stats.txt"),
            flagstat=path.join("results","samtools_qc","{sample}_samtools_flagstat.txt"),
            idxstats=path.join("results","samtools_qc","{sample}_samtools_idxstats.txt"),
        log: path.join("logs","samtools_qc","{sample}.log")
        threads: config["samtools_qc"].get("threads", 1)
        params: opt=config["samtools_qc"].get("opt", "")
        resources: mem_mb=config["samtools_qc"].get("mem", 1000)
        wrapper: "samtools_qc"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~COVERAGE RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "genomecov" in config:
    rule genomecov:
        input: rules.pbt_alignment_filter.output
        output: path.join("results","genomecov","{sample}.bedgraph")
        log: path.join("logs","genomecov","{sample}.log")
        threads: config["genomecov"].get("threads", 1)
        params: opt=config["genomecov"].get("opt", "")
        resources: mem_mb=config["genomecov"].get("mem", 1000)
        wrapper: "genomecov"

if "igvtools_count" in config:
    rule igvtools_count:
        input:
            bam = rules.pbt_alignment_filter.output,
            ref = config["reference"],
        output: path.join("results","igvtools_count","{sample}.tdf")
        log: path.join("logs","igvtools_count","{sample}.log")
        threads: config["igvtools_count"].get("threads", 1)
        params: opt=config["igvtools_count"].get("opt", "")
        resources: mem_mb=config["igvtools_count"].get("mem", 1000)
        wrapper: "igvtools_count"
