# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
output_tsv = "ref_CGI.tsv"
output_bed = "ref_CGI.bed"
output_bed_index = "ref_CGI.bed.idx"
output_tsv2 = "ref_CGI_2.tsv.gz"

# Rules
rule all:
    input: [output_tsv, output_bed, output_bed_index, output_tsv2]

rule pycometh_cgi_finder:
    input: ref=ref
    output: tsv=output_tsv, bed=output_bed, bed_index=output_bed_index
    log: "pycometh_cgi_finder.log"
    wrapper: "pycometh_cgi_finder"

rule pycometh_cgi_finder2:
    input: ref=ref
    output: tsv=output_tsv2
    params: opt="--merge_gap 50 --min_win_len 500 --min_CG_freq 0.45 --min_obs_CG_ratio 0.55"
    log: "pycometh_cgi_finder2.log"
    wrapper: "pycometh_cgi_finder"
