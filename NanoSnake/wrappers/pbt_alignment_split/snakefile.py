# Imports
from os.path import join

# Input and output data
bam_ont_input = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.bam")

chunks = [0,1,2,3]
bam_ont_output = expand("ont_reads_split_{c}.bam",c=chunks)
bam_illumina_output = expand("illumina_reads_split_{c}.bam",c=chunks)

# Rules
rule all:
    input: [bam_ont_output, bam_illumina_output]

rule alignment_split_ont:
    input: bam=bam_ont_input
    output: bam=bam_ont_output
    log: "alignment_split_ont.log"
    wrapper: "pbt_alignment_split"

rule alignment_split_illumina:
    input: bam=bam_illumina_input
    output: bam=bam_illumina_output
    log: "alignment_split_illumina.log"
    wrapper: "pbt_alignment_split"
