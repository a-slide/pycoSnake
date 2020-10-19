# Imports
from os.path import join

# Input and output data
input_vcf_1 = join(config["data_dir"], "ont_DNA", "sniffles_1.vcf")
input_vcf_2 = join(config["data_dir"], "ont_DNA", "sniffles_2.vcf")

output_vcf_1 = "filtered_reads_1.vcf"
output_vcf_2 = "filtered_reads_2.vcf"


# Rules
rule all:
    input: [output_vcf_1, output_vcf_2]

rule survivor_filter_1:
    input: vcf=input_vcf_1
    output: vcf=output_vcf_1
    threads: 4
    params: opt="50 500000 0.1 3"
    resources: mem_mb=1000
    log: "survivor_filter_1.log"
    wrapper: "survivor_filter"

rule survivor_filter_2:
    input: vcf=input_vcf_2
    output: vcf=output_vcf_2
    threads: 4
    params: opt="500 5000 0.2 5"
    resources: mem_mb=1000
    log: "survivor_filter_2.log"
    wrapper: "survivor_filter"
