# Imports
from os.path import join

# Input and output data
input_vcf_1 = join(config["data_dir"], "ont_DNA", "sniffles_1.vcf")
input_vcf_2 = join(config["data_dir"], "ont_DNA", "sniffles_2.vcf")
input_vcf_3 = join(config["data_dir"], "ont_DNA", "sniffles_3.vcf")

output_vcf = "merged_reads.vcf"

# Rules
rule all:
    input: [output_vcf]

rule survivor_merge:
    input: vcf = [input_vcf_1, input_vcf_2, input_vcf_3]
    output: vcf = output_vcf
    threads: 4
    params: opt="1000 1 1 -1 -1 -1"
    resources: mem_mb=1000
    log: "survivor_merge.log"
    wrapper: "survivor_merge"
