# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "small_ref.fa")
gff3 = join(config["data_dir"], "reference", "small_ref.gff3")
gtf = join(config["data_dir"], "reference", "small_ref.gtf")
index_gff3 = "index_gff3"
index_gtf = "index_gtf"

# Rules
rule all:
    input: [index_gff3, index_gtf]

rule star_index_gff3:
    input: ref=ref, annotation=gff3
    output: index_dir=directory(index_gff3)
    threads: 2
    params: opt="\
        --sjdbGTFfeatureExon exon\
        --sjdbGTFtagExonParentTranscript Parent\
        --sjdbGTFtagExonParentGene Name"
    resources: mem_mb = 1000
    log: "star_index_gff3.log"
    wrapper: "star_index"

rule star_index_gtf:
    input: ref=ref, annotation=gtf
    output: index_dir=directory(index_gtf)
    threads: 2
    params: opt=""
    resources: mem_mb=1000
    log: "star_index_gtf.log"
    wrapper: "star_index"
