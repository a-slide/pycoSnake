# Imports
from os.path import join

# Input and output data
bam_ont_input = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
bam_ont_output = "ont_reads.bam"
bam_illumina_output = "illumina_reads.bam"

# Rules
rule all:
    input: [bam_ont_output, bam_illumina_output]

rule alignment_filter_ont:
    input: bam=bam_ont_input
    output: bam=bam_ont_output, bam_index=bam_ont_output+".bai"
    log: "alignment_filter_ont.log"
    params: opt = "--min_align_len 1000 --min_freq_identity 0.7 --skip_unmapped --skip_secondary --skip_supplementary"
    wrapper: "pbt_alignment_filter"

rule alignment_filter_illumina:
    input: bam=bam_illumina_input
    output: bam=bam_illumina_output, bam_index=bam_illumina_output+".bai"
    log: "alignment_filter_illumina.log"
    params: opt = "--min_align_len 125 --min_freq_identity 0.8 --skip_unmapped --skip_secondary --skip_supplementary"
    wrapper: "pbt_alignment_filter"
