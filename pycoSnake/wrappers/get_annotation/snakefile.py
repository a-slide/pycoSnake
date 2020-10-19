# Imports
from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Input and output data
gff3_input = join(config["data_dir"], "reference", "ref.gff3.gz")
gff3_input_ftp = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-45/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.45.gff3.gz"
gff3_output_1 = "ref_1.gff3"
gtf_output_1 = "ref_1.gtf"
gff3_output_2 = "ref_2.gff3"
gtf_output_2 = "ref_2.gtf"

# Rules
rule all:
    input: [gff3_output_1, gtf_output_1, gff3_output_2, gtf_output_2]

rule get_annotation_from_local:
    input: gff3=gff3_input
    output: gff3=gff3_output_1, gtf=gtf_output_1
    params: opt=""
    log: "get_annotation_from_local.log"
    wrapper: "get_annotation"

rule get_annotation_from_ftp:
    input: gff3=FTP.remote(gff3_input_ftp)
    output: gff3=gff3_output_2, gtf=gtf_output_2
    params: opt=""
    log: "get_annotation_from_ftp.log"
    wrapper: "get_annotation"
