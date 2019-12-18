# -*- coding: utf-8 -*-

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Std lib
from os.path import join

# Third party lib
import pandas as pd

# Local imports
from NanoSnake.common import *
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~check config file version~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Minimum snakemake version
config_version=3
if not "config_version" in config or config["config_version"]!= config_version:
    raise NanoSnakeError ("Wrong configuration file version. Please regenerate config with -t config")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define samples sheet reference and getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list=list(sample_df.index)

def get_fastq (wildcards):
    return sample_df.loc[wildcards.sample, "fastq"]
def get_fast5 (wildcards):
    return sample_df.loc[wildcards.sample, "fast5"]
def get_seqsum (wildcards):
    return sample_df.loc[wildcards.sample, "seq_summary"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define IO for each rule~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build output files dictionnary
input_d=defaultdict(OrderedDict)
output_d=defaultdict(OrderedDict)
log_d=OrderedDict()

# Main rules
rule_name="preprocess_genone"
if config["genome"].startswith("ftp"):
    input_d[rule_name]["ref"]=FTP.remote(config["genome"])
elif config["genome"].startswith("http"):
    input_d[rule_name]["ref"]=HTTP.remote(config["genome"])
else:
    input_d[rule_name]["ref"]=config["genome"]
output_d[rule_name]["ref"]=join("results","main","genone","ref.fa")
output_d[rule_name]["index"]=join("results","main","genone","ref.fa.fai")
log_d[rule_name]=join("logs",rule_name,"ref.log")

rule_name="pbt_fastq_filter"
input_d[rule_name]["fastq"]=get_fastq
output_d[rule_name]["fastq"]=join("results","main","merged_fastq","{sample}.fastq")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="minimap2_index"
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
output_d[rule_name]["index"]=join("results","main","minimap2_index","ref.mmi")
log_d[rule_name]=join("logs",rule_name,"ref.log")

rule_name="minimap2_align"
input_d[rule_name]["index"]=output_d["minimap2_index"]["index"]
input_d[rule_name]["fastq"]=output_d["pbt_fastq_filter"]["fastq"]
output_d[rule_name]["bam"]=join("results","main","minimap2_alignments","{sample}.bam")
output_d[rule_name]["bam_index"]=join("results","main","minimap2_alignments","{sample}.bam.bai")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pbt_alignment_filter"
input_d[rule_name]["bam"]=output_d["minimap2_align"]["bam"]
output_d[rule_name]["bam"]=join("results","main","filtered_alignments","{sample}.bam")
output_d[rule_name]["bam_index"]=join("results","main","filtered_alignments","{sample}.bam.bai")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="nanopolish_index"
input_d[rule_name]["fastq"]=output_d["pbt_fastq_filter"]["fastq"]
input_d[rule_name]["fast5"]=get_fast5
input_d[rule_name]["seqsum"]=get_seqsum
output_d[rule_name]["index"]=join("results","main","merged_fastq","{sample}.fastq.index")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="nanopolish_call_methylation"
input_d[rule_name]["fastq"]=output_d["pbt_fastq_filter"]["fastq"]
input_d[rule_name]["fastq_index"]=output_d["nanopolish_index"]["index"]
input_d[rule_name]["bam"]=output_d["pbt_alignment_filter"]["bam"]
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","nanopolish_calls","{sample}.tsv")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycometh_cpg_aggregate"
input_d[rule_name]["tsv"]=output_d["nanopolish_call_methylation"]["tsv"]
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.tsv")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="ngmlr"
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
input_d[rule_name]["fastq"]=output_d["pbt_fastq_filter"]["fastq"]
output_d[rule_name]["bam"]=join("results","SV","ngmlr_alignments","{sample}.bam")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="sniffles"
input_d[rule_name]["bam"]=output_d["ngmlr"]["bam"]
output_d[rule_name]["vcf"]=join("results","SV","sniffles_calls","{sample}.vcf")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycoqc"
input_d[rule_name]["seqsum"]=get_seqsum
input_d[rule_name]["bam"]=output_d["minimap2_align"]["bam"]
output_d[rule_name]["html"]=join("results","QC","pycoqc","{sample}_pycoqc.html")
output_d[rule_name]["json"]=join("results","QC","pycoqc","{sample}_pycoqc.json")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="samtools_qc"
input_d[rule_name]["bam"]=output_d["minimap2_align"]["bam"]
output_d[rule_name]["stats"]=join("results","QC","samtools_qc","{sample}_samtools_stats.txt")
output_d[rule_name]["flagstat"]=join("results","QC","samtools_qc","{sample}_samtools_flagstat.txt")
output_d[rule_name]["idxstats"]=join("results","QC","samtools_qc","{sample}_samtools_idxstats.txt")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="bedtools_genomecov"
input_d[rule_name]["bam"]=output_d["pbt_alignment_filter"]["bam"]
output_d[rule_name]["bedgraph"]=join("results","coverage","bedgraph","{sample}.bedgraph")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="igvtools_count"
input_d[rule_name]["bam"]=output_d["pbt_alignment_filter"]["bam"]
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
output_d[rule_name]["tdf"]=join("results","coverage","igv_tdf","{sample}.tdf")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define all output depending on config file~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# main rules
run_rules = ["preprocess_genone", "pbt_fastq_filter", "minimap2_index", "minimap2_align", "pbt_alignment_filter"]
# optional methylation analysis
if all_in (config,["nanopolish_index", "nanopolish_call_methylation", "pycometh_cpg_aggregate"]):
    run_rules.extend (["nanopolish_index", "nanopolish_call_methylation", "pycometh_cpg_aggregate"])
# optional SV analysis
if all_in (config, ["ngmlr", "sniffles"]):
    run_rules.extend (["ngmlr", "sniffles"])
# Other optional orphan rules
for r in ["pycoqc", "samtools_qc", "bedtools_genomecov", "igvtools_count"]:
    if r in config:
        run_rules.append(r)

all_output=[]
for rule_name, rule_d in output_d.items():
    if rule_name in run_rules:
        for option, output in rule_d.items():
            all_output.append(output)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rule all:
    input: [expand(o, sample=sample_list) for o in all_output]

rule_name="preprocess_genone"
rule preprocess_genone:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "preprocess_genone"

rule_name="pbt_fastq_filter"
rule pbt_fastq_filter:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_fastq_filter"

rule_name="minimap2_index"
rule minimap2_index:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_index"

rule_name="minimap2_align"
rule minimap2_align:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name, 4)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_align"

rule_name="pbt_alignment_filter"
rule pbt_alignment_filter:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_filter"

rule_name="nanopolish_index"
rule nanopolish_index:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_index"

rule_name="nanopolish_call_methylation"
rule nanopolish_call_methylation:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name, 4)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_call_methylation"

rule_name="pycometh_cpg_aggregate"
rule pycometh_cpg_aggregate:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cpg_aggregate"

rule_name="ngmlr"
rule ngmlr:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name, 4)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "ngmlr"

rule_name="sniffles"
rule sniffles:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "sniffles"

rule_name="pycoqc"
rule pycoqc:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycoqc"

rule_name="samtools_qc"
rule samtools_qc:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "samtools_qc"

rule_name="bedtools_genomecov"
rule bedtools_genomecov:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "bedtools_genomecov"

rule_name="igvtools_count"
rule igvtools_count:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "igvtools_count"