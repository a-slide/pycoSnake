# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
annot = join(config["data_dir"], "reference", "cgi.bed")

input_tsv_1 = join(config["data_dir"], "ont_DNA", "Yeast_sample_1_CpG_aggregate.tsv.gz")
output_tsv_1 = "interval_aggregate_1.tsv"

input_tsv_2 = join(config["data_dir"], "ont_DNA", "Yeast_sample_2_CpG_aggregate.tsv.gz")
output_tsv_2 = "interval_aggregate_2.tsv.gz"
output_bed_2 = "interval_aggregate_2.bed"
output_bed_index_2 = "interval_aggregate_2.bed.idx"

# Rules
rule all:
    input: [output_tsv_1, output_tsv_2, output_bed_2, output_bed_index_2]

rule pycometh_interval_aggregate_1:
    input: tsv=input_tsv_1, ref=ref
    output: tsv=output_tsv_1
    log: "pycometh_interval_aggregate_1.log"
    params: opt="--interval_size 1000 --min_llr 1"
    wrapper: "pycometh_interval_aggregate"

rule pycometh_interval_aggregate_2:
    input: tsv=input_tsv_2, ref=ref, annot=annot
    output: tsv=output_tsv_2, bed=output_bed_2, bed_index=output_bed_index_2
    log: "pycometh_interval_aggregate_2.log"
    params:
        opt="--min_cpg_per_interval 3" ,
        sample_id="Test"
    wrapper: "pycometh_interval_aggregate"
