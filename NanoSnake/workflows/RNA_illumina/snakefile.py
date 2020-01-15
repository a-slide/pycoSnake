# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

def get_fastq1 (wildcards):
    return sample_df.loc[wildcards.sample, "fastq1"]
def get_fastq2 (wildcards):
    return sample_df.loc[wildcards.sample, "fastq2"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define IO for each rule~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build output files dictionnary
input_d=defaultdict(OrderedDict)
output_d=defaultdict(OrderedDict)
log_d=OrderedDict()

rule_name="preprocess_genone"
genome = config["genome"]
if genome.startswith("ftp"):
    input_d[rule_name]["ref"]=FTP.remote(genome)
elif genome.startswith("http"):
    input_d[rule_name]["ref"]=HTTP.remote(genome)
else:
    input_d[rule_name]["ref"]=genome
output_d[rule_name]["ref"]=join("results","main","genone","ref.fa")
output_d[rule_name]["index"]=join("results","main","genone","ref.fa.fai")
log_d[rule_name]=join("logs",rule_name,"preprocess_genone.log")

rule_name="preprocess_annotation"
annotation = config["annotation"]
if annotation.startswith("ftp"):
    input_d[rule_name]["annotation"]=FTP.remote(annotation)
elif annotation.startswith("http"):
    input_d[rule_name]["annotation"]=HTTP.remote(annotation)
else:
    input_d[rule_name]["annotation"]=annotation
output_d[rule_name]["gff3"]=join("results","main","annotation","ref.gff3")
output_d[rule_name]["gtf"]=join("results","main","annotation","ref.gtf")
log_d[rule_name]=join("logs",rule_name,"preprocess_annotation.log")

rule_name="fastp"
input_d[rule_name]["fastq1"]=get_fastq1
input_d[rule_name]["fastq2"]=get_fastq2
output_d[rule_name]["fastq1"]=join("results","main","fastp","{sample}_1.fastq.gz")
output_d[rule_name]["fastq2"]=join("results","main","fastp","{sample}_2.fastq.gz")
output_d[rule_name]["html"]=join("results","QC","fastp","{sample}_fastp.html")
output_d[rule_name]["json"]=join("results","QC","fastp","{sample}_fastp.json")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="star_index"
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
input_d[rule_name]["annotation"]=output_d["preprocess_annotation"]["gtf"]
output_d[rule_name]["index"]=join("results","main","star_index","SAindex")
log_d[rule_name]=join("logs",rule_name,"ref.log")

rule_name="star_align"
input_d[rule_name]["index"]=output_d["star_index"]["index"]
input_d[rule_name]["fastq1"]=output_d["fastp"]["fastq1"]
input_d[rule_name]["fastq2"]=output_d["fastp"]["fastq2"]
output_d[rule_name]["sj"]=join("results","main","star_alignments","{sample}_SJ.tsv")
output_d[rule_name]["bam"]=join("results","main","star_alignments","{sample}.bam")
output_d[rule_name]["count"]=join("results","counts","star","{sample}_counts.tsv")
output_d[rule_name]["bam_index"]=join("results","main","star_alignments","{sample}.bam.bai")
output_d[rule_name]["star_log"]=join("results","main","star_alignments","{sample}_Log.final.out")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pbt_alignment_filter"
input_d[rule_name]["bam"]=output_d["star_align"]["bam"]
output_d[rule_name]["bam"]=join("results","main","filtered_alignments","{sample}.bam")
output_d[rule_name]["bam_index"]=join("results","main","filtered_alignments","{sample}.bam.bai")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="star_count_merge"
input_d[rule_name]["counts"]=[join("results","counts","star",f"{sample}_counts.tsv") for sample in sample_list]
output_d[rule_name]["unstranded_counts"]=join("results","counts","star_merged","unstranded_counts.tsv")
output_d[rule_name]["positive_counts"]=join("results","counts","star_merged","positive_counts.tsv")
output_d[rule_name]["negative_counts"]=join("results","counts","star_merged","negative_counts.tsv")
log_d[rule_name]=join("logs",rule_name,"star_count_merge.log")

rule_name="cufflinks"
input_d[rule_name]["bam"]=output_d["star_align"]["bam"]
input_d[rule_name]["ref"]=output_d["preprocess_genone"]["ref"]
input_d[rule_name]["annotation"]=output_d["preprocess_annotation"]["gtf"]
output_d[rule_name]["genes_fpkm"]=join("results","counts","cufflinks","{sample}_genes_fpkm.tsv")
output_d[rule_name]["isoforms_fpkm"]=join("results","counts","cufflinks","{sample}_isoforms_fpkm.tsv")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="cufflinks_fpkm_merge"
input_d[rule_name]["fpkm_genes"]=[join("results","counts","cufflinks",f"{sample}_genes_fpkm.tsv") for sample in sample_list]
input_d[rule_name]["fpkm_isoforms"]=[join("results","counts","cufflinks",f"{sample}_isoforms_fpkm.tsv") for sample in sample_list]
output_d[rule_name]["fpkm_genes"]=join("results","counts","cufflinks_merged","fpkm_genes.tsv")
output_d[rule_name]["fpkm_isoforms"]=join("results","counts","cufflinks_merged","fpkm_isoforms.tsv")
log_d[rule_name]=join("logs",rule_name,"cufflinks_fpkm_merge.log")

rule_name="subread_featurecounts"
input_d[rule_name]["bam"]=output_d["star_align"]["bam"]
input_d[rule_name]["gtf"]=output_d["preprocess_annotation"]["gtf"]
output_d[rule_name]["counts"]=join("results","counts","featurecounts","{sample}_counts.tsv")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="subread_featurecounts_merge"
input_d[rule_name]["counts"]=[join("results","counts","featurecounts",f"{sample}_counts.tsv") for sample in sample_list]
output_d[rule_name]["counts"]=join("results","counts","featurecounts_merged","counts.tsv")
output_d[rule_name]["tpm"]=join("results","counts","featurecounts_merged","tpm.tsv")
log_d[rule_name]=join("logs",rule_name,"subread_featurecounts_merge.log")

rule_name="samtools_qc"
input_d[rule_name]["bam"]=output_d["star_align"]["bam"]
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
run_rules=["preprocess_genone", "preprocess_annotation", "fastp", "star_index", "star_align", "pbt_alignment_filter", "star_count_merge"]

# Optional cufflinks analysis
opt_rules = ["cufflinks", "cufflinks_fpkm_merge"]
if all_in(config, opt_rules):
    run_rules.extend(opt_rules)

# Optional subread_featurecounts analysis
opt_rules = ["subread_featurecounts", "subread_featurecounts_merge"]
if all_in(config, opt_rules):
    run_rules.extend(opt_rules)

# Other optional orphan rules
opt_rules = ["samtools_qc", "bedtools_genomecov", "igvtools_count"]
for r in opt_rules:
    if r in config:
        run_rules.extend(opt_rules)

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

rule_name="preprocess_annotation"
rule preprocess_annotation:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "preprocess_annotation"

rule_name="fastp"
rule fastp:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "fastp"

rule_name="star_index"
rule star_index:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "star_index"

rule_name="star_align"
rule star_align:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "star_align"

rule_name="pbt_alignment_filter"
rule pbt_alignment_filter:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_filter"

rule_name="star_count_merge"
rule star_count_merge:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "star_count_merge"

rule_name="cufflinks"
rule cufflinks:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "cufflinks"

rule_name="cufflinks_fpkm_merge"
rule cufflinks_fpkm_merge:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "cufflinks_fpkm_merge"

rule_name="subread_featurecounts"
rule subread_featurecounts:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "subread_featurecounts"

rule_name="subread_featurecounts_merge"
rule subread_featurecounts_merge:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "subread_featurecounts_merge"

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
