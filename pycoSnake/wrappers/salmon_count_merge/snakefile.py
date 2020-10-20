# Imports
from os.path import join
from glob import glob

# Input and output data
input_counts = sorted(glob(join(config["data_dir"], "illumina_RNA", "reads_*_salmon_quant.tsv")))
output_counts = "salmom_quant_counts.tsv"
output_tpm = "salmom_quant_tpm.tsv"

# Rules
rule all:
    input: [output_counts, output_tpm]

rule salmon_count_merge:
    input: counts=input_counts
    output: counts=output_counts, tpm=output_tpm
    log: "salmon_count_merge.log"
    wrapper: "salmon_count_merge"
