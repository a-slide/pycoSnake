# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Std lib
from os.path import join

# Third party lib
import pandas as pd
from snakemake.logging import logger

# Local imports
from pycoSnake.common import *
from snakemake.remote.FTP import RemoteProvider as FTP
from snakemake.remote.HTTP import RemoteProvider as HTTP

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def get_fastq (wildcards):
    return sample_df.loc[wildcards.sample, "fastq"]
def get_fast5 (wildcards):
    return sample_df.loc[wildcards.sample, "fast5"]
def get_seqsum (wildcards):
    return sample_df.loc[wildcards.sample, "seq_summary"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Initialise~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logger.info("Checking configuration file version")
config_version=12
if not "config_version" in config or config["config_version"]!= config_version:
    logger.error ("Wrong configuration file version. Please regenerate a template file with `--generate_template config -o`")
    sys.exit()
logger.debug(config)

logger.info("Loading sample file")
sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
if not sample_df.index.name=="sample_id" and not all_in(sample_df.columns, ["fastq","fast5","seq_summary"]):
    logger.error ("The provided sample sheet is in the correct format, Please regenerate a template file with `--generate_template sample_sheet -o`")
    sys.exit()
logger.debug(sample_df)
sample_list=list(sample_df.index)

logger.info("Define number of chunks")
try:
    nchunk=int(config["pbt_alignment_split"]["n_chunks"])
except:
    nchunk=4
chunk_list=list(range(nchunk))

logger.info("Specify way to download reference files")
ref=config["genome"]
if ref.startswith("ftp"):
    ref=FTP().remote(ref)
elif ref.startswith("http"):
    ref=HTTP().remote(ref)

gff3=config["annotation"]
if gff3.startswith("ftp"):
    gff3=FTP().remote(gff3)
elif gff3.startswith("http"):
    gff3=HTTP().remote(gff3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define all output depending on config file~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logger.info("Define conditional target files")
target_files=[]
target_files.extend(expand(join("results","main","filtered_alignments","{sample}.bam"), sample=sample_list))

if config["dna_methylation_call"] is True:
    logger.info("\tInclude rules for 'dna_methylation_call' module")
    target_files.extend(expand(join("results","methylation","nanopolish_calls","{sample}.tsv.gz"), sample=sample_list))

if config["differential_methylation"] is True:
    logger.info("\tInclude rules for 'differential_methylation' module")
    target_files.append(join("results","methylation","pycometh_comp_report", "pycoMeth_summary_report.html"))

if config["structural_variants_call"] is True:
    logger.info("\tInclude rules for 'structural_variants_call' module")
    target_files.append(join("results","SV","sniffles_all","merged.vcf"))

if config["quality_control"] is True:
    logger.info("\tInclude rules for 'quality_control' module")
    target_files.extend(expand(join("results","QC","pycoqc","{sample}_pycoqc.html"), sample=sample_list))
    target_files.extend(expand(join("results","QC","pycoqc","{sample}_pycoqc.json"), sample=sample_list))
    target_files.extend(expand(join("results","QC","samtools_qc","{sample}_samtools_stats.txt"), sample=sample_list))
    target_files.extend(expand(join("results","QC","samtools_qc","{sample}_samtools_flagstat.txt"), sample=sample_list))
    target_files.extend(expand(join("results","QC","samtools_qc","{sample}_samtools_idxstats.txt"), sample=sample_list))

if config["genome_coverage"] is True:
    logger.info("\tInclude rules for 'genome_coverage' module")
    target_files.extend(expand(join("results","coverage","bedgraph","{sample}.bedgraph"), sample=sample_list))
    target_files.extend(expand(join("results","coverage","igv_tdf","{sample}.tdf"), sample=sample_list))

for fn in target_files:
    logger.debug(f"\t{fn}")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input: target_files

rule_name="get_genome"
rule get_genome:
    input: ref=ref
    output:
        ref=join("results","input","genome","genome.fa"),
        index=join("results","input","genome","genome.fa.fai")
    log: join("logs",rule_name,"out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "get_genome"

rule_name="get_annotation"
rule get_annotation:
    input: gff3=gff3
    output:
        gff3=join("results","input","annotation","annotation.gff3"),
        gtf=join("results","input","annotation","annotation.gtf")
    log: join("logs", rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "get_annotation"

rule_name="pbt_fastq_filter"
rule pbt_fastq_filter:
    input: fastq=get_fastq
    output: fastq=join("results","input","merged_fastq","{sample}.fastq")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_fastq_filter"

rule_name="minimap2_index"
rule minimap2_index:
    input: ref=rules.get_genome.output.ref
    output: index=join("results","main","minimap2_index","ref.mmi")
    log: join("logs",rule_name,"ref.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_index"

rule_name="minimap2_align"
rule minimap2_align:
    input:
        index=rules.minimap2_index.output.index,
        fastq=rules.pbt_fastq_filter.output.fastq
    output:
        bam=join("results","main","minimap2_alignments","{sample}.bam"),
        bam_index=join("results","main","minimap2_alignments","{sample}.bam.bai")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "minimap2_align"

rule_name="pbt_alignment_filter"
rule pbt_alignment_filter:
    input: bam=rules.minimap2_align.output.bam
    output:
        bam=join("results","main","filtered_alignments","{sample}.bam"),
        bam_index=join("results","main","filtered_alignments","{sample}.bam.bai")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_filter"

rule_name="nanopolish_index"
rule nanopolish_index:
    input:
        fastq=rules.pbt_fastq_filter.output.fastq,
        fast5=get_fast5,
        seqsum=get_seqsum
    output: index=join("results","input","merged_fastq","{sample}.fastq.index")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_index"

rule_name="pbt_alignment_split"
rule pbt_alignment_split:
    input: bam=rules.minimap2_align.output.bam
    output:
        bam=temp(expand(join("results","methylation","split_alignments","{{sample}}","{chunk}.bam"), chunk=chunk_list)),
        bam_index=temp(expand(join("results","methylation","split_alignments","{{sample}}","{chunk}.bam.bai"), chunk=chunk_list))
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pbt_alignment_split"

rule_name="nanopolish_call_methylation"
rule nanopolish_call_methylation:
    input:
        fastq=rules.pbt_fastq_filter.output.fastq,
        fastq_index=rules.nanopolish_index.output.index,
        bam=join("results","methylation","split_alignments","{sample}","{chunk}.bam"),
        bam_index=join("results","methylation","split_alignments","{sample}","{chunk}.bam.bai"),
        ref=rules.get_genome.output.ref
    output: tsv=temp(join("results","methylation","nanopolish_calls","{sample}","{chunk}.tsv"))
    log: join("logs",rule_name,"{sample}","{chunk}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_call_methylation"

rule_name="nanopolish_concat"
rule nanopolish_concat:
    input: tsv_list=expand(join("results","methylation","nanopolish_calls","{{sample}}","{chunk}.tsv"), chunk=chunk_list)
    output: tsv=protected(join("results","methylation","nanopolish_calls","{sample}.tsv.gz"))
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "nanopolish_concat"

rule_name="pycometh_cgi_finder"
rule pycometh_cgi_finder:
    input: ref=rules.get_genome.output.ref
    output:
        tsv=join("results","methylation","pycometh_cgi_finder","CGI.tsv.gz"),
        bed=join("results","methylation","pycometh_cgi_finder","CGI.bed"),
        bed_index=join("results","methylation","pycometh_cgi_finder","CGI.bed.idx")
    log: join("logs",rule_name,"ref.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cgi_finder"

rule_name="pycometh_cpg_aggregate"
rule pycometh_cpg_aggregate:
    input:
        tsv= rules.nanopolish_concat.output.tsv,
        ref=rules.get_genome.output.ref
    output:
        tsv=join("results","methylation","pycometh_cpg_aggregate","{sample}.tsv.gz"),
        bed=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed"),
        bed_index=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed.idx")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cpg_aggregate"

rule_name="pycometh_interval_aggregate"
rule pycometh_interval_aggregate:
    input:
        tsv=rules.pycometh_cpg_aggregate.output.tsv,
        ref=rules.get_genome.output.ref,
        annot=rules.pycometh_cgi_finder.output.bed
    output:
        tsv=join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz"),
        bed=join("results","methylation","pycometh_interval_aggregate","{sample}.bed"),
        bed_index=join("results","methylation","pycometh_interval_aggregate","{sample}.bed.idx")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_interval_aggregate"

rule_name="pycometh_meth_comp"
rule pycometh_meth_comp:
    input:
        tsv=expand(join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz"), sample=sample_list),
        ref=rules.get_genome.output.ref
    output:
        tsv=join("results","methylation","pycometh_meth_comp","meth_comp.tsv.gz"),
        bed=join("results","methylation","pycometh_meth_comp","meth_comp.bed"),
        bed_index=join("results","methylation","pycometh_meth_comp","meth_comp.bed.idx")
    log: join("logs",rule_name,"meth_comp.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_meth_comp"

rule_name="pycometh_comp_report"
rule pycometh_comp_report:
    input:
        tsv=rules.pycometh_meth_comp.output.tsv,
        gff3=rules.get_annotation.output.gff3,
        ref=rules.get_genome.output.ref
    output: summary_report=join("results","methylation","pycometh_comp_report", "pycoMeth_summary_report.html")
    log: join("logs",rule_name,"comp_report.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_comp_report"

rule_name="ngmlr"
rule ngmlr:
    input:
        ref=rules.get_genome.output.ref,
        fastq=rules.pbt_fastq_filter.output.fastq
    output: bam=join("results","SV","ngmlr_alignments","{sample}.bam")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "ngmlr"

rule_name="sniffles"
rule sniffles:
    input: bam=rules.ngmlr.output.bam
    output: vcf=join("results","SV","sniffles","{sample}_raw.vcf")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "sniffles"

rule_name="survivor_filter"
rule survivor_filter:
    input: vcf=rules.sniffles.output.vcf
    output: vcf=join("results","SV","sniffles","{sample}_filtered.vcf")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "survivor_filter"

rule_name="survivor_merge"
rule survivor_merge:
    input: vcf=expand(join("results","SV","sniffles","{sample}_filtered.vcf"), sample=sample_list)
    output: vcf=join("results","SV","sniffles","merged.vcf")
    log: join("logs",rule_name,"merged.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "survivor_merge"

rule_name="sniffles_all"
rule sniffles_all:
    input:
        bam=rules.ngmlr.output.bam,
        vcf=rules.survivor_merge.output.vcf
    output: vcf=join("results","SV","sniffles_all","{sample}_raw.vcf")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "sniffles"

rule_name="survivor_merge_all"
rule survivor_merge_all:
    input: vcf=expand(join("results","SV","sniffles_all","{sample}_raw.vcf"), sample=sample_list)
    output: vcf=join("results","SV","sniffles_all","merged.vcf")
    log: join("logs",rule_name,"merged.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "survivor_merge"

rule_name="pycoqc"
rule pycoqc:
    input:
        seqsum=get_seqsum,
        bam=rules.minimap2_align.output.bam
    output:
        html=join("results","QC","pycoqc","{sample}_pycoqc.html"),
        json=join("results","QC","pycoqc","{sample}_pycoqc.json")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycoqc"

rule_name="samtools_qc"
rule samtools_qc:
    input: bam=rules.minimap2_align.output.bam
    output:
        stats=join("results","QC","samtools_qc","{sample}_samtools_stats.txt"),
        flagstat=join("results","QC","samtools_qc","{sample}_samtools_flagstat.txt"),
        idxstats=join("results","QC","samtools_qc","{sample}_samtools_idxstats.txt")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "samtools_qc"

rule_name="bedtools_genomecov"
rule bedtools_genomecov:
    input: bam=rules.pbt_alignment_filter.output.bam
    output: bedgraph=join("results","coverage","bedgraph","{sample}.bedgraph")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "bedtools_genomecov"

rule_name="igvtools_count"
rule igvtools_count:
    input:
        bam=rules.pbt_alignment_filter.output.bam,
        ref=rules.get_genome.output.ref
    output: tdf=join("results","coverage","igv_tdf","{sample}.tdf")
    log: join("logs",rule_name,"{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "igvtools_count"
