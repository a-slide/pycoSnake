# Imports
from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Input and output data
ref_input = "ftp://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.28.cdna.all.fa.gz"
ref_output = "ref.fa"
index_dir = "index"

# Rules
rule all:
    input: [ref_output, index_dir]

rule get_transcriptome:
    input: ref=FTP.remote(ref_input)
    output: ref=ref_output
    log: "get_transcriptome.log"
    wrapper: "get_transcriptome"

rule salmon_index:
    input: ref=ref_output
    output: index_dir=directory(index_dir)
    params: opt=""
    threads: 4
    log: "salmon_index.log"
    wrapper: "salmon_index"
