# Imports
from os.path import join

# Input and output data
bam_ngmlr = join(config["data_dir"], "ont_DNA", "reads_ngmlr.bam")
vcf = "reads.vcf"

# Rules
rule all:
    input: [vcf]

rule sniffles_1:
    input: bam=bam_ngmlr
    output: vcf=vcf
    threads: 4
    params: opt="--min_support 1 --max_num_splits 7 --min_length 30"
    resources: mem_mb=1000
    log: "sniffles_1.log"
    wrapper: "sniffles"
