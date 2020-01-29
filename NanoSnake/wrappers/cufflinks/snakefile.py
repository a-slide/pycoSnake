# Imports
from os.path import join

# Input and output data
bam_illumina_1 = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
bam_illumina_2 = join(config["data_dir"], "illumina_RNA", "reads_2.bam")
ref = join(config["data_dir"], "reference", "small_ref.fa")
gff3 = join(config["data_dir"], "reference", "small_ref.gff3")
gtf = join(config["data_dir"], "reference", "small_ref.gtf")
genes_fpkm_1 = "genes_fpkm_1.tsv"
isoforms_fpkm_1 = "isoforms_fpkm_1.tsv"
genes_fpkm_2 = "genes_fpkm_2.tsv"
isoforms_fpkm_2 = "isoforms_fpkm_2.tsv"
transcript_gtf_2 = "transcript_gtf_2.gtf"

# Rules
rule all:
    input: [genes_fpkm_1, isoforms_fpkm_1, genes_fpkm_2, transcript_gtf_2]

rule cufflinks_1_gff3:
    input: bam=bam_illumina_1, ref=ref, annotation=gff3
    output: genes_fpkm=genes_fpkm_1, isoforms_fpkm=isoforms_fpkm_1
    log: "cufflinks_1.log"
    params: opt = "--library-type fr-firststrand --upper-quartile-norm"
    threads: 4
    resources: mem_mb = 1000
    wrapper: "cufflinks"

rule cufflinks_2_gtf:
    input: bam=bam_illumina_2, ref=ref, annotation=gtf
    output: genes_fpkm=genes_fpkm_2, transcript_gtf=transcript_gtf_2, isoforms_fpkm=isoforms_fpkm_2
    log: "cufflinks_2.log"
    params: opt = "--library-type fr-firststrand"
    threads: 4
    resources: mem_mb = 1000
    wrapper: "cufflinks"
