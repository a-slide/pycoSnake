# Imports
from os.path import join
from glob import glob

# Input and output data
input_counts=counts = sorted(glob(join(config["data_dir"], "illumina_RNA", "reads_*_featurecounts.tsv")))
output_counts_1 = "featurecounts_merge_counts_1.tsv"
output_tpm_2 = "featurecounts_merge_tpm_2.tsv"
output_counts_3 = "featurecounts_merge_counts_3.tsv"
output_tpm_3 = "featurecounts_merge_tpm_3.tsv"

# Rules
rule all:
    input: [output_counts_1, output_tpm_2, output_counts_3, output_tpm_3]

rule subread_featurecounts_merge_1:
    input: counts=input_counts
    output: counts=output_counts_1
    log: "subread_featurecounts_merge_1.log"
    wrapper: "subread_featurecounts_merge"

rule subread_featurecounts_merge_2:
    input: counts=input_counts
    output: tpm=output_tpm_2
    log: "subread_featurecounts_merge_2.log"
    wrapper: "subread_featurecounts_merge"

rule subread_featurecounts_merge_3:
    input: counts=input_counts
    output: counts=output_counts_3, tpm=output_tpm_3
    log: "subread_featurecounts_merge_3.log"
    wrapper: "subread_featurecounts_merge"
