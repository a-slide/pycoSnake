# Imports
from os.path import join
from glob import glob

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
input_tsv_1 = sorted(glob(join(config["data_dir"], "ont_DNA", "Yeast_sample_*_CpG_aggregate.tsv.gz")))
input_tsv_2 = sorted(glob(join(config["data_dir"], "ont_DNA", "Yeast_sample_*_Interval_aggregate.tsv.gz")))
output_bed_1 = "meth_comp_cpg.bed"
output_bed_idx_1 = "meth_comp_cpg.bed.idx"
output_tsv_1 = "meth_comp_cpg.tsv.gz"
output_bed_2 = "meth_comp_interval_1.bed"
output_tsv_2 = "meth_comp_interval_1.tsv.gz"
output_bed_3 = "meth_comp_interval_2.bed"
output_tsv_3 = "meth_comp_interval_2.tsv.gz"

# Rules
rule all:
    input: [output_bed_1, output_tsv_1, output_bed_idx_1, output_bed_2, output_tsv_2, output_bed_3, output_tsv_3]

rule pycometh_meth_comp_1:
    input: tsv=input_tsv_1, ref=ref
    output: bed=output_bed_1, tsv=output_tsv_1, bed_index=output_bed_idx_1
    log: "pycometh_meth_comp_1.log"
    wrapper: "pycometh_meth_comp"

rule pycometh_meth_comp_2:
    input: tsv=input_tsv_2, ref=ref
    output: bed=output_bed_2, tsv=output_tsv_2
    log: "pycometh_meth_comp_2.log"
    wrapper: "pycometh_meth_comp"

rule pycometh_meth_comp_3:
    input: tsv=input_tsv_2, ref=ref
    output: bed=output_bed_3, tsv=output_tsv_3
    params: opt="--max_missing 2 --min_diff_llr 0 --pvalue_threshold 0.1"
    log: "pycometh_meth_comp_3.log"
    wrapper: "pycometh_meth_comp"
