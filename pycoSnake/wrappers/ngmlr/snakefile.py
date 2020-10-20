# Imports
from os.path import join

# Input and output data
ref_input = join(config["data_dir"], "reference", "ref.fa.gz")
fastq = join(config["data_dir"], "ont_DNA", "reads.fastq.gz")
fasta = join(config["data_dir"], "ont_DNA", "sim_read_1.fa.gz")

ref_output = "ref.fa"
bam_1 = "reads_ngmlr_1.bam"
bam_index_1 = bam_1+".bai"
bam_2 = "reads_ngmlr_2.bam"
bam_index_2 = bam_2+".bai"

# Rules
rule all:
    input: [ref_output, bam_1, bam_index_1, bam_2, bam_index_2]

# Copy genome because ngmlr add an index on the fly
rule get_genone:
    input: ref=ref_input
    output: ref=ref_output
    log: "get_genone.log"
    wrapper: "get_genome"

rule ngmlr_from_fq:
    input: fastq=fastq, ref=ref_output
    output: bam=bam_1, bam_index=bam_index_1
    threads: 4
    params: opt="-x ont"
    resources: mem_mb=1000
    log: "ngmlr_from_fq.log"
    wrapper: "ngmlr"

rule ngmlr_from_fa:
    input: fastq=fasta, ref=ref_output
    output: bam=bam_2, bam_index=bam_index_2
    threads: 4
    params: opt="-x ont"
    resources: mem_mb=1000
    log: "ngmlr_from_fa.log"
    wrapper: "ngmlr"
