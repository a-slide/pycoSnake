# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "small_ref.fa")
fastq = join(config["data_dir"], "ont_DNA", "reads.fastq.gz")
index = "ref.mmi"
bam = "reads.bam"
bam_index = "reads.bam.bai"

# Rules
rule all:
    input: [index, bam, bam_index]

rule minimap2_index:
    input: ref=ref
    output: index=index
    threads: 2
    params: opt=""
    resources: mem_mb=1000
    log: "minimap2_index.log"
    wrapper: "minimap2_index"

rule minimap2_align:
    input: fastq=fastq, index=index
    output: bam=bam, bam_index=bam_index
    threads: 2
    params: opt="-x map-ont -L"
    resources: mem_mb=1000
    log: "minimap2_align.log"
    wrapper: "minimap2_align"
