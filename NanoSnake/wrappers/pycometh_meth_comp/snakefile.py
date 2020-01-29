# Imports
from os.path import join
from glob import glob

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
input_tsv_1 = sorted(glob(join(config["data_dir"], "ont_DNA", "Yeast_sample_*_CpG_aggregate.tsv.gz")))
output_bed_1 = "meth_comp_1.bed"
output_bed_index_1 = "meth_comp_1.bed.idx"
output_tsv_1 = "meth_comp_2.tsv.gz"

# Rules
rule all:
    input: [output_bed_1, output_tsv_1]

rule pycometh_meth_comp_1:
    input: tsv=input_tsv_1, ref=ref
    output: bed=output_bed_1, tsv=output_tsv_1, bed_index=output_bed_index_1
    params: opt="--max_missing 1 --min_diff_llr 1"
    log: "pycometh_meth_comp_1.log"
    wrapper: "pycometh_meth_comp"
