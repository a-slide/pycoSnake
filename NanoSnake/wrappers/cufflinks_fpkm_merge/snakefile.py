# Imports
from os.path import join
from glob import glob

# Input and output data
input_gene_fpkm = sorted(glob(join(config["data_dir"], "illumina_RNA", "genes_fpkm_*.tsv")))
input_isoform_fpkm = sorted(glob(join(config["data_dir"], "illumina_RNA", "isoforms_fpkm_*.tsv")))
output_gene_fpkm_1 = "cufflinks_gene_fpkm_1.tsv"
output_isoform_fpkm_1 = "cufflinks_isoform_fpkm_1.tsv"
output_isoform_fpkm_2 = "cufflinks_isoform_fpkm_2.tsv"

# Rules
rule all:
    input: [output_gene_fpkm_1, output_isoform_fpkm_1, output_isoform_fpkm_2]

rule cufflinks_fpkm_merge_1:
    input: fpkm_genes=input_gene_fpkm, fpkm_isoforms=input_isoform_fpkm
    output: fpkm_genes=output_gene_fpkm_1, fpkm_isoforms=output_isoform_fpkm_1
    log: "cufflinks_fpkm_merge_1.log"
    wrapper: "cufflinks_fpkm_merge"

rule cufflinks_fpkm_merge_2:
    input: fpkm_isoforms=input_isoform_fpkm
    output: fpkm_isoforms=output_isoform_fpkm_2
    log: "cufflinks_fpkm_merge_2.log"
    wrapper: "cufflinks_fpkm_merge"
