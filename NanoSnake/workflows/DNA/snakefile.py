# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Std lib
from os.path import join as pj

# Third party lib
import pandas as pd
from snakemake.utils import min_version

# Local imports
from NanoSnake.common import *
from NanoSnake import __version__ as package_version

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Preliminary checks and imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Minimum snakemake version
min_version("5.4.2")

# Minimum snakemake version
config_version = 2
if not "config_version" in config or config["config_version"]!= config_version:
    raise NanoSnakeError ("Wrong configuration file version. Please regenerate config with -g")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define samples sheet and getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sample_df = pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list = list(sample_df.index)

def get_fastq (wildcards):
    return sample_df.loc[wildcards.sample, "fastq"]
def get_fast5 (wildcards):
    return sample_df.loc[wildcards.sample, "fast5"]
def get_seq_summary (wildcards):
    return sample_df.loc[wildcards.sample, "seq_summary"]

#~~~~~~~~~~~~~~~~~~~~~~~~Define directories subworflow flags and shortcuts~~~~~~~~~~~~~~~~~~~~~~~~~#

res_dir ={
    "pbt_fastq_filter": pj("results","main","merged_fastq"),
    "minimap2_index": pj("results","main","minimap2_index"),
    "minimap2_align": pj("results","main","minimap2_alignments"),
    "pbt_alignment_filter": pj("results","main","filtered_alignments"),
    "nanopolish_index": pj("results","main","merged_fastq"),
    "nanopolish_call_methylation": pj("results","methylation","nanopolish_calls"),
    "pycometh_aggregate": pj("results","methylation","pycometh_aggregate"),
    "ngmlr": pj("results","SV","ngmlr_alignments"),
    "sniffles": pj("results","SV","sniffles_calls"),
    "pycoqc": pj("results","QC","pycoqc"),
    "samtools_qc": pj("results","QC","samtools_qc"),
    "genomecov": pj("results","coverage","bedgraph"),
    "igvtools_count": pj("results","coverage","igv")}

log_dir ={
    "pbt_fastq_filter": pj("logs","pbt_fastq_filter"),
    "minimap2_index": pj("logs","minimap2_index"),
    "minimap2_align": pj("logs","minimap2_align"),
    "pbt_alignment_filter": pj("logs","pbt_alignment_filter"),
    "nanopolish_index": pj("logs","nanopolish_index"),
    "nanopolish_call_methylation": pj("logs","nanopolish_call_methylation"),
    "pycometh_aggregate": pj("logs","pycometh_aggregate"),
    "ngmlr": pj("logs","ngmlr"),
    "sniffles": pj("logs","sniffles"),
    "pycoqc": pj("logs","pycoqc"),
    "samtools_qc": pj("logs","samtools_qc"),
    "genomecov": pj("logs","genomecov"),
    "igvtools_count": pj("logs","igvtools_count")}

run_meth_calling = all_in (config, ["nanopolish_index", "nanopolish_call_methylation", "pycometh_aggregate"])
run_SV_calling = all_in (config, ["ngmlr", "sniffles"])
run_pycoqc = "pycoqc" in config
run_samtools_qc = "samtools_qc" in config
run_genomecov = "genomecov" in config
run_igvtools_count = "igvtools_count" in config
reference = config["reference"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Top Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input:
        expand(pj(res_dir["pbt_fastq_filter"],"{sample}.fastq"), sample=sample_list),
        expand(pj(res_dir["minimap2_index"],"ref.mmi")),
        expand(pj(res_dir["minimap2_align"], "{sample}.bam"), sample=sample_list),
        expand(pj(res_dir["pbt_alignment_filter"], "{sample}.bam"), sample=sample_list),
        expand(pj(res_dir["pbt_fastq_filter"],"{sample}.fastq.index"), sample=sample_list) if run_meth_calling else [],
        expand(pj(res_dir["nanopolish_call_methylation"],"{sample}.tsv"), sample=sample_list) if run_meth_calling else [],
        expand(pj(res_dir["pycometh_aggregate"],"{sample}.bed"), sample=sample_list) if run_meth_calling else [],
        expand(pj(res_dir["pycometh_aggregate"],"{sample}.tsv"), sample=sample_list) if run_meth_calling else [],
        expand(pj(res_dir["ngmlr"],"{sample}.bam"), sample=sample_list) if run_SV_calling else [],
        expand(pj(res_dir["sniffles"],"{sample}.vcf"), sample=sample_list) if run_SV_calling else [],
        expand(pj(res_dir["pycoqc"],"{sample}_pycoqc.html"), sample=sample_list) if run_pycoqc else [],
        expand(pj(res_dir["pycoqc"],"{sample}_pycoqc.json"), sample=sample_list) if run_pycoqc else [],
        expand(pj(res_dir["samtools_qc"],"{sample}_samtools_stats.txt"), sample=sample_list) if run_samtools_qc else [],
        expand(pj(res_dir["samtools_qc"],"{sample}_samtools_flagstat.txt"), sample=sample_list) if run_samtools_qc else [],
        expand(pj(res_dir["samtools_qc"],"{sample}_samtools_idxstats.txt"), sample=sample_list) if run_samtools_qc else [],
        expand(pj(res_dir["genomecov"],"{sample}.bedgraph"), sample=sample_list) if run_genomecov else [],
        expand(pj(res_dir["igvtools_count"],"{sample}.tdf"), sample=sample_list) if run_igvtools_count else [],

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Main rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule_name = "pbt_fastq_filter"
rule pbt_fastq_filter:
    input: get_fastq
    output: pj(res_dir[rule_name],"{sample}.fastq")
    log: pj(log_dir[rule_name],"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_fastq_filter"

rule_name = "minimap2_index"
rule minimap2_index:
    input: reference
    output: pj(res_dir[rule_name],"ref.mmi")
    log: pj(log_dir[rule_name],"ref.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_index"

rule_name = "minimap2_align"
rule minimap2_align:
    input:
        index=rules.minimap2_index.output,
        fastq=rules.pbt_fastq_filter.output
    output: pj(res_dir[rule_name],"{sample}.bam")
    log: pj(log_dir[rule_name],"{sample}.log")
    threads: get_threads(config, rule_name, 4)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_align"

rule_name = "pbt_alignment_filter"
rule pbt_alignment_filter:
    input: rules.minimap2_align.output
    output: pj(res_dir[rule_name],"{sample}.bam")
    log: pj(log_dir[rule_name],"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_filter"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~dna methylation rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if run_meth_calling:

    rule_name = "nanopolish_index"
    rule nanopolish_index:
        input:
            fastq = rules.pbt_fastq_filter.output,
            fast5 = get_fast5,
            seq_summary = get_seq_summary,
        output: pj(res_dir[rule_name],"{sample}.fastq.index")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name),
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "nanopolish_index"

    rule_name = "nanopolish_call_methylation"
    rule nanopolish_call_methylation:
        input:
            fastq = rules.pbt_fastq_filter.output,
            fastq_index = rules.nanopolish_index.output,
            bam = rules.pbt_alignment_filter.output,
            ref = reference,
        output: pj(res_dir[rule_name],"{sample}.tsv")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name, 4)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "nanopolish_call_methylation"

    rule_name = "pycometh_aggregate"
    rule pycometh_aggregate:
        input:
            call = rules.nanopolish_call_methylation.output,
            ref = reference,
        output:
            bed = pj(res_dir[rule_name],"{sample}.bed"),
            tsv = pj(res_dir[rule_name],"{sample}.tsv"),
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "pycometh_aggregate"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SV rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if run_SV_calling:

    rule_name = "ngmlr"
    rule ngmlr:
        input:
            ref=reference,
            fastq=rules.pbt_fastq_filter.output
        output: pj(res_dir[rule_name],"{sample}.bam")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name, 4)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "ngmlr"

    rule_name = "sniffles"
    rule sniffles:
        input: rules.ngmlr.output
        output: pj(res_dir[rule_name],"{sample}.vcf")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "sniffles"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QC rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if run_pycoqc:
    rule_name = "pycoqc"
    rule pycoqc:
        input:
            seq_summary = get_seq_summary,
            bam = rules.minimap2_align.output
        output:
            html=pj(res_dir[rule_name],"{sample}_pycoqc.html"),
            json=pj(res_dir[rule_name],"{sample}_pycoqc.json")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "pycoqc"

if run_samtools_qc:
    rule_name = "samtools_qc"
    rule samtools_qc:
        input: rules.minimap2_align.output,
        output:
            stats=pj(res_dir[rule_name],"{sample}_samtools_stats.txt"),
            flagstat=pj(res_dir[rule_name],"{sample}_samtools_flagstat.txt"),
            idxstats=pj(res_dir[rule_name],"{sample}_samtools_idxstats.txt"),
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "samtools_qc"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~coverage rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if run_genomecov:
    rule_name = "genomecov"
    rule genomecov:
        input: rules.pbt_alignment_filter.output
        output: pj(res_dir[rule_name],"{sample}.bedgraph")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "genomecov"

if run_igvtools_count:
    rule_name = "igvtools_count"
    rule igvtools_count:
        input:
            bam = rules.pbt_alignment_filter.output,
            ref = reference,
        output: pj(res_dir[rule_name],"{sample}.tdf")
        log: pj(log_dir[rule_name],"{sample}.log")
        threads: get_threads(config, rule_name)
        params: opt=get_opt(config, rule_name)
        resources: mem_mb=get_mem(config, rule_name)
        wrapper: "igvtools_count"
