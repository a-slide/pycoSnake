# Imports
from os.path import join

# Input and output data
fast5 = join(config["data_dir"], "ont_DNA", "fast5")
input_fastq = join(config["data_dir"], "ont_DNA", "reads.fastq.gz")
seqsum = join(config["data_dir"], "ont_DNA", "sequencing_summary.txt")
output_fastq = "reads.fastq"
index = output_fastq+".index"

# Rules
rule all:
    input: [output_fastq, index]

# Copy and extract fastq
rule pbt_fastq_filter:
    input: fastq=input_fastq
    output: fastq=output_fastq
    log: "pbt_fastq_filter.log"
    wrapper: "pbt_fastq_filter"

# Index reads
rule nanopolish_index:
    input: fast5=fast5, fastq=output_fastq, seqsum=seqsum
    output: index=index
    threads: 4
    params: opt=""
    resources: mem_mb=1000
    log: "nanopolish_index.log"
    wrapper: "nanopolish_index"
