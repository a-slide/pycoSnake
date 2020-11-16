# Imports
from os.path import join
from glob import glob

# Input and output data
gff3 = join(config["data_dir"], "reference", "ref.gff3")
ref = join(config["data_dir"], "reference", "ref.fa")
input_tsv = join(config["data_dir"], "ont_DNA", "Yeast_meth_comp.tsv.gz")
summary_report_1 = join("Comp_Report_Yeast_1", "pycoMeth_summary_report.html")
summary_report_2 = join("Comp_Report_Yeast_2", "pycoMeth_summary_report.html")

# Rules
rule all:
    input: [summary_report_1, summary_report_2]

rule pycometh_comp_report_1:
    input: tsv=input_tsv, gff3=gff3, ref=ref
    output: summary_report=summary_report_1
    log: "pycometh_comp_report_1.log"
    wrapper: "pycometh_comp_report"

rule pycometh_comp_report_2:
    input: tsv=input_tsv, gff3=gff3, ref=ref
    output: summary_report=summary_report_2
    params: opt="--max_tss_distance 10000 --pvalue_threshold 1 --min_diff_llr 0"
    log: "pycometh_comp_report_2.log"
    wrapper: "pycometh_comp_report"
