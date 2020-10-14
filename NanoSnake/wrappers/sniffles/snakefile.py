# Imports
from os.path import join

# Input and output data
bam_1 = join(config["data_dir"], "ont_DNA", "ngmlr_sim_read_1.bam")
bam_2 = join(config["data_dir"], "ont_DNA", "ngmlr_sim_read_2.bam")
bam_3 = join(config["data_dir"], "ont_DNA", "ngmlr_sim_read_3.bam")
input_vcf = join(config["data_dir"], "ont_DNA", "sniffles_merged.vcf")

output_vcf_1 = "sniffles_1.vcf"
output_vcf_2 = "sniffles_2.vcf"
output_vcf_3 = "sniffles_3.vcf"

# Rules
rule all:
    input: [output_vcf_1, output_vcf_2, output_vcf_3]

rule sniffles_1:
    input: bam=bam_1
    output: vcf=output_vcf_1
    threads: 2
    params: opt="--min_support 3 --max_num_splits 7 --max_distance 1000 --min_length 100 --minmapping_qual 20 --min_seq_size 1000 --allelefreq 0.1"
    resources: mem_mb=1000
    log: "sniffles_1.log"
    wrapper: "sniffles"

rule sniffles_2:
    input: bam=bam_2
    output: vcf=output_vcf_2
    threads: 2
    params: opt=""
    resources: mem_mb=1000
    log: "sniffles_2.log"
    wrapper: "sniffles"

rule sniffles_3:
    input: bam=bam_3, vcf=input_vcf
    output: vcf=output_vcf_3
    threads: 2
    params: opt="--min_support 5 --max_num_splits 7 --max_distance 1000 --min_length 50 --minmapping_qual 20 --min_seq_size 500 --allelefreq 0.2"
    resources: mem_mb=1000
    log: "sniffles_3.log"
    wrapper: "sniffles"
