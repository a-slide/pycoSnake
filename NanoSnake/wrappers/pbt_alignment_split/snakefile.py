# Imports
from os.path import join

# Input and output data
bam_ont_input = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.bam")

chunks = [0,1,2,3]
bam_ont_output = [f"ont_reads_split_{{sample}}_{c}.bam" for c in chunks]
bam_illumina_output = "illumina_reads_split_{c}.bam"

# Rules
rule all:
    input: [bam_ont_output, expand(bam_illumina_output, c=chunks)]

rule alignment_split_ont:
    input: bam=bam_ont_input
    output: bam=bam_ont_output
    log: "alignment_split_ont.log"
    wrapper: "pbt_alignment_split"

checkpoint alignment_split_illumina:
    input: bam=bam_illumina_input
    output: bam=directory
    log: "alignment_split_illumina.log"
    wrapper: "pbt_alignment_split"
