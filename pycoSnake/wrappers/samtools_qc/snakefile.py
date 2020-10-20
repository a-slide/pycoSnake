# Imports
from os.path import join

# Input and output data
bam_ont = join(config["data_dir"], "ont_DNA", "reads.bam")
bam_illumina = join(config["data_dir"], "illumina_RNA", "reads_1.bam")
stats_ont = "stats_ont.txt"
flagstat_ont = "flagstat_ont.txt"
idxstats_ont = "idxstats_ont.txt"
stats_illumina = "stats_illumina.txt"
flagstat_illumina = "flagstat_illumina.txt"
idxstats_illumina = "idxstats_illumina.txt"

# Rules
rule all:
    input: [stats_ont, flagstat_ont, idxstats_ont, stats_illumina, flagstat_illumina, idxstats_illumina]

rule samtools_qc_ont:
    input: bam=bam_ont
    output: stats=stats_ont ,flagstat=flagstat_ont ,idxstats=idxstats_ont
    log: "samtools_qc_ont.log"
    wrapper: "samtools_qc"

rule samtools_qc_illumina:
    input: bam=bam_illumina
    output: stats=stats_illumina ,flagstat=flagstat_illumina ,idxstats=idxstats_illumina
    log: "samtools_qc_illumina.log"
    wrapper: "samtools_qc"
