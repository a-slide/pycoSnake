# Imports
from os.path import join

# Input and output data
bam_1 = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
bam_2 = join(config["data_dir"], "illumina_RNA", "reads_2.bam")
bam_3 = join(config["data_dir"], "ont_DNA", "reads.bam")
gtf = join(config["data_dir"], "reference", "small_ref.gtf")
counts_1 = "counts_1.tsv"
counts_2 = "counts_2.tsv"
counts_3 = "counts_3.tsv"

# Rules
rule all:
    input: [counts_1, counts_2, counts_3]

rule featurecounts_1:
    input: bam=bam_1, gtf=gtf
    output: counts=counts_1
    params: opt = "-p"
    threads: 4
    resources: mem_mb = 1000
    log: "featurecounts_1.log"
    wrapper: "subread_featurecounts"

rule featurecounts_2:
    input: bam=bam_2, gtf=gtf
    output: counts=counts_2
    params: opt = "-p"
    threads: 4
    resources: mem_mb = 1000
    log: "featurecounts_2.log"
    wrapper: "subread_featurecounts"

rule featurecounts_3_ont:
    input: bam=bam_3, gtf=gtf
    output: counts=counts_3
    params: opt = "-L"
    threads: 1
    resources: mem_mb = 1000
    log: "featurecounts_3.log"
    wrapper: "subread_featurecounts"
