# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
input_tsv_1 = join(config["data_dir"], "ont_DNA", "nanopolish_sample_1.tsv.gz")
output_tsv_1 = "cpg_aggregate_1.tsv"
input_tsv_2 = join(config["data_dir"], "ont_DNA", "nanopolish_sample_2.tsv.gz")
output_tsv_2 = "cpg_aggregate_2.tsv.gz"
output_bed_2 = "cpg_aggregate_2.bed"
output_bed_index_2 = "cpg_aggregate_2.bed.idx"

# Rules
rule all:
    input: [output_tsv_1, output_tsv_2, output_bed_2, output_bed_index_2]

rule pycometh_cpg_aggregate_1:
    input: tsv=input_tsv_1, ref=ref
    output: tsv=output_tsv_1
    log: "pycometh_cpg_aggregate_1.log"
    wrapper: "pycometh_cpg_aggregate"

rule pycometh_cpg_aggregate_2:
    input: tsv=input_tsv_2, ref=ref
    output: tsv=output_tsv_2, bed=output_bed_2, bed_index=output_bed_index_2
    log: "pycometh_cpg_aggregate_2.log"
    params:
        opt="--min_depth 5 --min_llr 1",
        sample_id="Test"
    wrapper: "pycometh_cpg_aggregate"
