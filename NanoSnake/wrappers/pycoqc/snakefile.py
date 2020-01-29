# Imports
from os.path import join

# Input and output data
seqsum = join(config["data_dir"], "ont_DNA", "sequencing_summary.txt")
bam = join(config["data_dir"], "ont_DNA", "reads.bam")
html_1 = "pycoqc_report_1.html"
json_1 = "pycoqc_report_1.json"
html_2 = "pycoqc_report_2.html"
json_2 = "pycoqc_report_2.json"

# Rules
rule all:
    input: [html_1, json_1, html_2, json_2]

rule pycoqc_noBAM:
    input: seqsum=seqsum
    output: html=html_1, json=json_1
    params: opt = "--filter_calibration --filter_duplicated --min_pass_len 100 --min_pass_qual 7", sample_id = "S1"
    log: "pycoqc_1.log"
    wrapper: "pycoqc"

rule pycoqc_BAM:
    input: seqsum=seqsum, bam=bam
    output: html=html_2, json=json_2
    params: opt = "--filter_calibration --filter_duplicated --min_pass_len 100 --min_pass_qual 7", sample_id = "S2"
    log: "pycoqc_2.log"
    wrapper: "pycoqc"
