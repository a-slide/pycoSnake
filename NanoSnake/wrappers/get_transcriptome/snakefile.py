# Imports
from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Input and output data
ref1_input = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
ref1_output = "ref1.fa"
ref1_tsv = "ref1.tsv"
ref1_index = ref1_output+".fai"
ref2_input = "ftp://ftp.ensembl.org/pub/release-99/fasta/oryzias_latipes/cdna/Oryzias_latipes.ASM223467v1.cdna.all.fa.gz"
ref2_output = "ref2.fa"
ref2_tsv = "ref2.tsv"
ref2_index = ref2_output+".fai"
ref3_input = "ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.abinitio.fa.gz"
ref3_output = "ref3.fa"
ref3_tsv = "ref3.tsv"
ref3_index = ref3_output+".fai"

# Rules
rule all:
    input: [ref1_output, ref1_index, ref1_tsv, ref2_output, ref2_index, ref2_tsv, ref3_output, ref3_index, ref3_tsv]

rule get_transcriptome_yeast:
    input: ref=FTP.remote(ref1_input)
    output: ref=ref1_output, tsv=ref1_tsv, index=ref1_index
    log: "get_transcriptome_yeast.log"
    wrapper: "get_transcriptome"

rule get_transcriptome_medaka:
    input: ref=FTP.remote(ref2_input)
    output: ref=ref2_output, tsv=ref2_tsv, index=ref2_index
    log: "get_transcriptome_medaka.log"
    wrapper: "get_transcriptome"

rule get_transcriptome_mus:
    input: ref=FTP.remote(ref3_input)
    output: ref=ref3_output, tsv=ref3_tsv, index=ref3_index
    log: "get_transcriptome_mus.log"
    wrapper: "get_transcriptome"
