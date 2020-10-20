# Imports
from os.path import join

# Input and output data
fastq_ont_input = join(config["data_dir"], "ont_DNA", "reads.fastq.gz")
fastq_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.fastq.gz")
fastq_ont_output = "fastq_ont_output.gz"
fastq_illumina_output = "fastq_illumina_output.gz"

# Rules
rule all:
    input: [fastq_ont_output, fastq_illumina_output]

rule fastq_filter_ont:
    input: fastq=fastq_ont_input
    output: fastq=fastq_ont_output
    log: "fastq_filter_ont.log"
    params: opt = "--remove_duplicates --min_len 100 --min_qual 7"
    wrapper: "pbt_fastq_filter"

rule fastq_filter_illumina:
    input: fastq=fastq_illumina_input
    output: fastq=fastq_illumina_output
    log: "fastq_filter_illumina.log"
    params: opt = "--remove_duplicates --min_qual 20"
    wrapper: "pbt_fastq_filter"
