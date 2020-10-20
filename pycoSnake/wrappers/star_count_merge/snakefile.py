# Imports
from os.path import join
from glob import glob

# Input and output data
counts = sorted(glob(join(config["data_dir"], "illumina_RNA", "reads_*_counts.tsv")))
unstranded_counts_1 = "star_unstranded_counts_1.tsv"
positive_counts_1 = "star_positive_counts_1.tsv"
unstranded_counts_2 = "star_unstranded_counts_2.tsv"
negative_counts_2 = "star_negative_counts_1.tsv"

# Rules
rule all:
    input: [unstranded_counts_1, positive_counts_1, unstranded_counts_2, negative_counts_2]

rule star_count_merge_1:
    input: counts=counts
    output: unstranded_counts=unstranded_counts_1, positive_counts=positive_counts_1
    log: "star_count_merge_1.log"
    wrapper: "star_count_merge"

rule star_count_merge_2:
    input: counts=counts
    output: unstranded_counts=unstranded_counts_2, negative_counts=negative_counts_2
    log: "star_count_merge_2.log"
    wrapper: "star_count_merge"
