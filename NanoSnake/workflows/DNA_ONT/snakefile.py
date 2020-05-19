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
config_version=8
if not "config_version" in config or config["config_version"]!= config_version:
    raise NanoSnakeError ("Wrong configuration file version. Please regenerate config with `--generate_template config -o`")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define samples sheet reference and getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list=list(sample_df.index)

def get_fastq (wildcards):
    return sample_df.loc[wildcards.sample, "fastq"]
def get_fast5 (wildcards):
    return sample_df.loc[wildcards.sample, "fast5"]
def get_seqsum (wildcards):
    return sample_df.loc[wildcards.sample, "seq_summary"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define number of chunks for nanopolish~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
try:
    nchunk = int(config["pbt_alignment_split"]["n_chunks"])
except:
    nchunk = 2
chunk_list = list(range(nchunk))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define IO for each rule~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build output files dictionnary
input_d=defaultdict(OrderedDict)
output_d=defaultdict(OrderedDict)
log_d=OrderedDict()

rule_name="get_genome"
ref = config["genome"]
if ref.startswith("ftp"): input_d[rule_name]["ref"]=FTP.remote(ref)
elif ref.startswith("http"): input_d[rule_name]["ref"]=HTTP.remote(ref)
else: input_d[rule_name]["ref"]=ref
output_d[rule_name]["ref"]=join("results","input","genome","genome.fa")
output_d[rule_name]["index"]=join("results","input","genome","genome.fa.fai")
log_d[rule_name]=join("logs",rule_name,"get_genome.log")

rule_name="get_annotation"
gff3 = config["annotation"]
if gff3.startswith("ftp"): input_d[rule_name]["gff3"]=FTP.remote(gff3)
elif gff3.startswith("http"): input_d[rule_name]["gff3"]=HTTP.remote(gff3)
else: input_d[rule_name]["gff3"]=gff3
output_d[rule_name]["gff3"]=join("results","input","annotation","annotation.gff3")
output_d[rule_name]["gtf"]=join("results","input","annotation","annotation.gtf")
log_d[rule_name]=join("logs",rule_name,"get_annotation.log")

rule_name="pbt_fastq_filter"
input_d[rule_name]["fastq"]=get_fastq
output_d[rule_name]["fastq"]=join("results","input","merged_fastq","{sample}.fastq")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="minimap2_index"
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
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
output_d[rule_name]["index"]=join("results","input","merged_fastq","{sample}.fastq.index")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pbt_alignment_split"
input_d[rule_name]["bam"]=output_d["pbt_alignment_filter"]["bam"]
output_d[rule_name]["bam"]=expand(join("results","methylation","split_alignments","{{sample}}","{chunk}.bam"), chunk=chunk_list) # Expand chunks
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="nanopolish_call_methylation"
input_d[rule_name]["fastq"]=output_d["pbt_fastq_filter"]["fastq"]
input_d[rule_name]["fastq_index"]=output_d["nanopolish_index"]["index"]
input_d[rule_name]["bam"]=join("results","methylation","split_alignments","{sample}","{chunk}.bam")
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","nanopolish_calls","{sample}","{chunk}.tsv")
log_d[rule_name]=join("logs",rule_name,"{sample}","{chunk}.log")

rule_name="pycometh_cgi_finder"
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_cgi_finder","CGI.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_cgi_finder","CGI.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_cgi_finder","CGI.bed.idx")
log_d[rule_name]=join("logs",rule_name,"ref.log")

rule_name="pycometh_cpg_aggregate"
input_d[rule_name]["tsv"]=expand(join("results","methylation","nanopolish_calls","{{sample}}","{chunk}.tsv"), chunk=chunk_list) # Aggregate chunks
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed.idx")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycometh_interval_aggregate"
input_d[rule_name]["tsv"]=output_d["pycometh_cpg_aggregate"]["tsv"]
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
input_d[rule_name]["annot"]=output_d["pycometh_cgi_finder"]["bed"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_interval_aggregate","{sample}.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_interval_aggregate","{sample}.bed.idx")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycometh_meth_comp"
input_d[rule_name]["tsv"]=expand(join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz"), sample=sample_list) # Aggregate samples
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_meth_comp","meth_comp.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_meth_comp","meth_comp.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_meth_comp","meth_comp.bed.idx")
log_d[rule_name]=join("logs",rule_name,"meth_comp.log")

rule_name="pycometh_comp_report"
input_d[rule_name]["tsv"]=output_d["pycometh_meth_comp"]["tsv"]
input_d[rule_name]["gff3"]=output_d["get_annotation"]["gff3"]
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["outdir"]=directory(join("results","methylation","pycometh_comp_report"))
log_d[rule_name]=join("logs",rule_name,"comp_report.log")

rule_name="ngmlr"
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
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
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tdf"]=join("results","coverage","igv_tdf","{sample}.tdf")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define all output depending on config file~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# main rules
run_rules = ["get_genome", "get_annotation", "pbt_fastq_filter", "minimap2_index", "minimap2_align", "pbt_alignment_filter"]

# Optional nanopolish analysis analysis
np_rules = ["pbt_alignment_split", "nanopolish_index", "nanopolish_call_methylation"]
if all_in (config, np_rules):
    run_rules.extend (np_rules)

# Optional pycoMeth analysis
pcm_rules = ["pycometh_cgi_finder", "pycometh_cpg_aggregate", "pycometh_interval_aggregate", "pycometh_meth_comp", "pycometh_comp_report"]
if all_in (config, np_rules) and all_in (config, pcm_rules):
    run_rules.extend (pcm_rules)

# Optional SV analysis
SV_rules = ["ngmlr", "sniffles"]
if all_in (config, SV_rules):
    run_rules.extend (SV_rules)

# Other optional orphan rules
orphan_rules = ["pycoqc", "samtools_qc", "bedtools_genomecov", "igvtools_count"]
for r in orphan_rules:
    if r in config:
        run_rules.append(r)

all_output=[]
for rule_name, rule_d in output_d.items():
    if rule_name in run_rules:
        for option, output in rule_d.items():
            all_output.append(output)
all_output = flatten_list(all_output)

all_expand = []
for output in all_output:
    all_expand.append(list(set(expand(output, sample=sample_list, chunk=chunk_list))))
all_expand = flatten_list(all_expand)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input: all_expand

rule_name="get_genome"
rule get_genome:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "get_genome"

rule_name="get_annotation"
rule get_annotation:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "get_annotation"

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

rule_name="pbt_alignment_split"
rule pbt_alignment_split:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_split"

rule_name="nanopolish_call_methylation"
rule nanopolish_call_methylation:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name, 4)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_call_methylation"

rule_name="pycometh_cgi_finder"
rule pycometh_cgi_finder:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cgi_finder"

rule_name="pycometh_cpg_aggregate"
rule pycometh_cpg_aggregate:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cpg_aggregate"

rule_name="pycometh_interval_aggregate"
rule pycometh_interval_aggregate:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_interval_aggregate"

rule_name="pycometh_meth_comp"
rule pycometh_meth_comp:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_meth_comp"

rule_name="pycometh_comp_report"
rule pycometh_comp_report:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_comp_report"

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
