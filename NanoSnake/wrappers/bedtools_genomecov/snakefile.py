# Imports
from os.path import join

# Input and output data
bam_ont_input = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
bg_ont_output = "ont_reads.bedgraph"
bg_illumina_output = "illumina_reads.bedgraph"

# Rules
rule all:
    input: [bg_ont_output, bg_illumina_output]

rule genomecov_ont:
    input: bam=bam_ont_input
    output: bedgraph=bg_ont_output
    log: "genomecov_ont.log"
    params: opt="", sample_id="ont"
    wrapper: "bedtools_genomecov"

rule genomecov_illumina:
    input: bam=bam_illumina_input
    output: bedgraph=bg_illumina_output
    log: "genomecov_illumina.log"
    params: opt="", sample_id="illumina"
    wrapper: "bedtools_genomecov"
