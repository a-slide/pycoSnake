# Imports
from os.path import join

# Input and output data
input_gz1 = join(config["data_dir"], "ont_DNA", "nanopolish_sample_1.tsv.gz")
input_gz2 = join(config["data_dir"], "ont_DNA", "nanopolish_sample_2.tsv.gz")
input_nc1 = join(config["data_dir"], "ont_DNA", "sniffles_1.vcf")
input_nc2 = join(config["data_dir"], "ont_DNA", "sniffles_1.vcf")

output_1 = "merged_1.txt"
output_2 = "merged_2.txt.gz"
output_3 = "merged_3.txt.gz"

# Rules
rule all:
    input: [output_1, output_2, output_3]

rule nanopolish_concat_1:
    input: tsv_list=[input_gz1, input_gz2]
    output: tsv=output_1
    log: "nanopolish_concat_1.log"
    wrapper: "nanopolish_concat"

rule nanopolish_concat_2:
    input: tsv_list=[input_nc1, input_nc2]
    output: tsv=output_2
    log: "nanopolish_concat_2.log"
    wrapper: "nanopolish_concat"

rule nanopolish_concat_3:
    input: tsv_list=[input_gz1, input_gz2, input_nc1, input_nc2]
    output: tsv=output_3
    log: "nanopolish_concat_3.log"
    wrapper: "nanopolish_concat"
