# Imports
from os.path import join

# Input and output data
ref_input = join(config["data_dir"], "reference", "ref.fa.gz")
fastq = join(config["data_dir"], "ont_DNA", "reads.fastq.gz")
ref_output = "ref.fa"
bam = "reads_ngmlr.bam"
bam_index = bam+".bai"

# Rules
rule all:
    input: [ref_output, bam, bam_index]

# Copy genome because ngmlr add an index on the fly
rule get_genone:
    input: ref=ref_input
    output: ref=ref_output
    log: "preprocess_genone.log"
    wrapper: "get_genome"

rule ngmlr:
    input: fastq=fastq, ref=ref_output
    output: bam=bam, bam_index=bam_index
    threads: 4
    params: opt="-x ont"
    resources: mem_mb=1000
    log: "ngmlr.log"
    wrapper: "ngmlr"
