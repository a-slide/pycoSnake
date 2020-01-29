# Imports
from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Input and output data
ref_input = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
fastq1 = join(config["data_dir"], "illumina_RNA", "reads_1.fastq.gz")
fastq2 = join(config["data_dir"], "illumina_RNA", "reads_2.fastq.gz")
fastq3 = join(config["data_dir"], "illumina_RNA", "reads_3.fastq.gz")
fastq4 = join(config["data_dir"], "illumina_RNA", "reads_4.fastq.gz")
ref_output = "ref.fa"
index_dir = "index"
quant_dir_1 = "quant_res_1"
counts_1 = "salmon_quant_1.tsv"
quant_dir_2 = "quant_res_2"
counts_2 = "salmon_quant_2.tsv"

# Rules
rule all:
    input: [ref_output, index_dir, quant_dir_1, counts_1, quant_dir_2, counts_2]

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

rule salmon_quant_1:
    input: index_dir=index_dir, fastq1=fastq1, fastq2=fastq2
    output: quant_dir=directory(quant_dir_1), counts=counts_1
    params: opt="--libType A --validateMappings"
    threads: 4
    log: "salmon_quant_1.log"
    wrapper: "salmon_quant"

rule salmon_quant_2:
    input: index_dir=index_dir, fastq1=fastq3, fastq2=fastq4
    output: quant_dir=directory(quant_dir_2), counts=counts_2
    params: opt="--libType IU"
    threads: 4
    log: "salmon_quant_2.log"
    wrapper: "salmon_quant"
