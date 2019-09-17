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
        expand(path.join("results","merge_filter_fastq","{sample}.fastq"), sample=sample_list),
        expand(path.join("results","minimap2_index","ref.mmi")),
        expand(path.join("results","minimap2_align", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","samtools_filter", "{sample}.bam"), sample=sample_list),
        expand(path.join("results","merge_filter_fastq","{sample}.fastq.index"), sample=sample_list) if "nanopolish_index" in config else [],
        expand(path.join("results", "nanopolish_call_methylation","{sample}.tsv"), sample=sample_list) if "nanopolish_call_methylation" in config else [],
        expand(path.join("results", "nanopolishcomp_freq_meth_calculate","{sample}.bed"), sample=sample_list) if "nanopolishcomp_freq_meth_calculate" in config else [],
        expand(path.join("results", "nanopolishcomp_freq_meth_calculate","{sample}.tsv"), sample=sample_list) if "nanopolishcomp_freq_meth_calculate" in config else [],
        expand(path.join("results","pycoqc","{sample}_pycoqc.html"), sample=sample_list) if "pycoqc" in config else [],
        expand(path.join("results","pycoqc","{sample}_pycoqc.json"), sample=sample_list) if "pycoqc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_stats.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_flagstat.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","samtools_qc","{sample}_samtools_idxstats.txt"), sample=sample_list) if "samtools_qc" in config else [],
        expand(path.join("results","genomecov","{sample}.bedgraph"), sample=sample_list) if "genomecov" in config else [],
        expand(path.join("results","igvtools_count","{sample}.tdf"), sample=sample_list) if "igvtools_count" in config else [],

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CORE RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule merge_filter_fastq:
    input: get_fastq
    output: path.join("results","merge_filter_fastq","{sample}.fastq")
    log: path.join("logs","merge_filter_fastq","{sample}.log")
    threads: config["merge_filter_fastq"].get("threads", 1)
    params: opt=config["merge_filter_fastq"].get("opt", "")
    resources: mem_mb=config["merge_filter_fastq"].get("mem", 1000)
    wrapper: "merge_filter_fastq"

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

rule samtools_filter:
    input: rules.minimap2_align.output
    output: path.join("results","samtools_filter","{sample}.bam")
    log: path.join("logs","samtools_filter","{sample}.log")
    threads: config["samtools_filter"].get("threads", 2)
    params: opt=config["samtools_filter"].get("opt", "")
    resources: mem_mb=config["samtools_filter"].get("mem", 1000)
    wrapper: "samtools_filter"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DNA METHYLATION RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "nanopolish_index" in config:
    rule nanopolish_index:
        input:
            fastq = rules.merge_filter_fastq.output,
            fast5 = get_fast5,
            seq_summary = get_seq_summary,
        output: path.join("results","merge_filter_fastq","{sample}.fastq.index")
        log: path.join("logs","nanopolish_index","{sample}.log")
        threads: config["nanopolish_index"].get("threads", 2)
        params: opt=config["nanopolish_index"].get("opt", ""),
        resources: mem_mb=config["nanopolish_index"].get("mem", 1000)
        wrapper: "nanopolish_index"

if "nanopolish_call_methylation" in config:
    rule nanopolish_call_methylation:
        input:
            fastq = rules.merge_filter_fastq.output,
            fastq_index = rules.nanopolish_index.output,
            bam = rules.samtools_filter.output,
            ref = config["reference"],
        output: path.join("results", "nanopolish_call_methylation","{sample}.tsv")
        log: path.join("logs","nanopolish_call_methylation","{sample}.log")
        threads: config["nanopolish_call_methylation"].get("threads", 2)
        params: opt=config["nanopolish_call_methylation"].get("opt", ""),
        resources: mem_mb=config["nanopolish_call_methylation"].get("mem", 1000)
        wrapper: "nanopolish_call_methylation"

if "nanopolishcomp_freq_meth_calculate" in config:
    rule nanopolishcomp_freq_meth_calculate:
        input:
            call = rules.nanopolish_call_methylation.output,
            ref = config["reference"],
        output:
            bed = path.join("results", "nanopolishcomp_freq_meth_calculate","{sample}.bed"),
            tsv = path.join("results", "nanopolishcomp_freq_meth_calculate","{sample}.tsv"),
        log: path.join("logs","nanopolishcomp_freq_meth_calculate","{sample}.log")
        threads: config["nanopolishcomp_freq_meth_calculate"].get("threads", 1)
        params: opt=config["nanopolishcomp_freq_meth_calculate"].get("opt", ""),
        resources: mem_mb=config["nanopolishcomp_freq_meth_calculate"].get("mem", 1000)
        wrapper: "nanopolishcomp_freq_meth_calculate"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QC RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if "pycoqc" in config:
    rule pycoqc:
        input:
            seq_summary = get_seq_summary
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
        input: rules.samtools_filter.output
        output: path.join("results","genomecov","{sample}.bedgraph")
        log: path.join("logs","genomecov","{sample}.log")
        threads: config["genomecov"].get("threads", 1)
        params: opt=config["genomecov"].get("opt", "")
        resources: mem_mb=config["genomecov"].get("mem", 1000)
        wrapper: "genomecov"

if "igvtools_count" in config:
    rule igvtools_count:
        input:
            bam = rules.samtools_filter.output,
            ref = config["reference"],
        output: path.join("results","igvtools_count","{sample}.tdf")
        log: path.join("logs","igvtools_count","{sample}.log")
        threads: config["igvtools_count"].get("threads", 1)
        params: opt=config["igvtools_count"].get("opt", "")
        resources: mem_mb=config["igvtools_count"].get("mem", 1000)
        wrapper: "igvtools_count"
