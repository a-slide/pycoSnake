# Imports
from os.path import join

# Input and output data
fastq1_in = join(config["data_dir"], "illumina_RNA", "reads_1.fastq.gz")
fastq2_in = join(config["data_dir"], "illumina_RNA", "reads_2.fastq.gz")
fastq1_pe_out = "illumina_reads_pe_1.fastq.gz"
fastq2_pe_out = "illumina_reads_pe_2.fastq.gz"
json_pe = "illumina_reads_pe.json"
html_pe = "illumina_reads_pe.html"
fastq1_se_out = "illumina_reads_se_1.fastq.gz"
json_se = "illumina_reads_se.json"
html_se = "illumina_reads_se.html"

# Rules
rule all:
    input: [fastq1_pe_out, fastq2_pe_out, html_pe, json_pe, fastq1_se_out, json_se, html_se]

rule fastp_pe:
    input: fastq1=fastq1_in, fastq2=fastq2_in
    output: fastq1=fastq1_pe_out, fastq2=fastq2_pe_out, html=html_pe, json=json_pe
    log: "fastp_pe.log"
    threads: 2
    params: opt = ""
    resources: mem_mb = 1000
    wrapper: "fastp"

rule fastp_se:
    input: fastq1=fastq1_in
    output: fastq1=fastq1_se_out, json=json_se, html=html_se
    log: "fastp_se.log"
    threads: 2
    params: opt = ""
    resources: mem_mb = 1000
    wrapper: "fastp"
