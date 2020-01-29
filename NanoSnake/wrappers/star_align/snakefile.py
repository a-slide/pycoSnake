# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
gff3 = join(config["data_dir"], "reference", "ref.gff3")
fastq1 = join(config["data_dir"], "illumina_RNA", "reads_1.fastq.gz")
fastq2 = join(config["data_dir"], "illumina_RNA", "reads_2.fastq.gz")
fastq3 = join(config["data_dir"], "illumina_RNA", "reads_3.fastq.gz")
fastq4 = join(config["data_dir"], "illumina_RNA", "reads_4.fastq.gz")
index_dir = "index"
sj_1="illumina_reads_1_SJ.tsv"
count_1="illumina_reads_1_counts.tsv"
bam_1="illumina_reads_1.bam"
count_2="illumina_reads_2_counts.tsv"
star_log_2="illumina_reads_2_star.log"


# Rules
rule all:
    input: [sj_1, count_1, bam_1, count_2, star_log_2]

rule star_index:
    input: ref=ref, annotation=gff3
    output: index_dir=directory(index_dir)
    threads: 2
    params: opt=" \
        --sjdbGTFfeatureExon exon\
        --sjdbGTFtagExonParentTranscript Parent\
        --sjdbGTFtagExonParentGene Name"
    resources: mem_mb = 1000
    log: "star_index.log"
    wrapper: "star_index"

rule star_align_1:
    input: index_dir=index_dir, fastq1=fastq1, fastq2=fastq2
    output: sj=sj_1, count=count_1, bam=bam_1, bam_index=bam_1+".bai"
    threads: 4
    params: opt = ""
    resources: mem_mb = 1000
    log: "star_align_1.log"
    wrapper: "star_align"

rule star_align_2:
    input: index_dir=index_dir, fastq1=fastq3, fastq2=fastq4
    output: count=count_2, star_log=star_log_2
    threads: 4
    params: opt = "\
        --outFilterType BySJout\
        --outFilterMultimapNmax 20\
        --alignSJoverhangMin 8\
        --alignSJDBoverhangMin 1\
        --outFilterMismatchNmax 999\
        --outFilterMismatchNoverLmax 0.04\
        --alignIntronMin 20\
        --alignIntronMax 1000000\
        --alignMatesGapMax 1000000"
    resources: mem_mb = 1000
    log: "star_align_2.log"
    wrapper: "star_align"
