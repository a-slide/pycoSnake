# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
bam_ont_input = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina_input = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
tdf_ont_output = "ont_reads.tdf"
tdf_illumina_output = "illumina_reads.tdf"

# Rules
rule all:
    input: [tdf_ont_output, tdf_illumina_output]

rule igvtools_count_ont:
    input: bam=bam_ont_input, ref=ref
    output: tdf=tdf_ont_output
    log: "igvtools_count_ont.log"
    params: opt=""
    wrapper: "igvtools_count"

rule igvtools_count_illumina:
    input: bam=bam_illumina_input, ref=ref
    output: tdf=tdf_illumina_output
    log: "igvtools_count_illumina.log"
    params: opt=""
    wrapper: "igvtools_count"
