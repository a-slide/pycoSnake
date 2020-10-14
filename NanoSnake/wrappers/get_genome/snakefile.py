# Imports
from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

# Input and output data
ref1_input = join(config["data_dir"], "reference", "small_ref.fa")
ref2_input = join(config["data_dir"], "reference", "ref.fa.gz")
ref3_input = "ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.MT.fa.gz"
ref1_output = "ref1.fa"
ref1_index = ref1_output+".fai"
ref2_output = "ref2.fa"
ref2_index = ref2_output+".fai"
ref3_output = "ref3.fa"
ref3_index = ref3_output+".fai"

# Rules
rule all:
    input: [ref1_output, ref1_index, ref2_output, ref2_index, ref3_output, ref3_index]

rule get_genome_from_fa:
    input: ref=ref1_input
    output: ref=ref1_output, index=ref1_index
    log: "get_genome_from_fa.log"
    wrapper: "get_genome"

rule get_genome_from_gz:
    input: ref=ref2_input
    output: ref=ref2_output, index=ref2_index
    log: "get_genome_from_gz.log"
    wrapper: "get_genome"

rule get_genome_from_ftp:
    input: ref=FTP.remote(ref3_input)
    output: ref=ref3_output, index=ref3_index
    log: "get_genome_from_ftp.log"
    wrapper: "get_genome"
