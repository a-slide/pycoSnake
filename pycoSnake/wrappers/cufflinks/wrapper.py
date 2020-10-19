__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.2"

# Imports
from snakemake.shell import shell
import tempfile
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
bam = snakemake.input.bam
ref = snakemake.input.ref
annotation = snakemake.input.annotation
genes_fpkm = snakemake.output.get("genes_fpkm", None)
isoforms_fpkm = snakemake.output.get("isoforms_fpkm", None)
transcript_gtf = snakemake.output.get("transcript_gtf", None)
outdir = os.path.dirname(os.path.abspath(snakemake.output[0]))

# Run shell command
with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:
    shell("cufflinks {opt} -p {snakemake.threads} -G {annotation} -b {ref} -o {temp_dir} {bam} &> {snakemake.log}")

    # Only keep fpkm files
    if genes_fpkm:
        temp_genes_fpkm = os.path.join(temp_dir, "genes.fpkm_tracking")
        shell("mv {temp_genes_fpkm} {genes_fpkm}")
    if isoforms_fpkm:
        temp_isoforms_fpkm = os.path.join(temp_dir, "isoforms.fpkm_tracking")
        shell("mv {temp_isoforms_fpkm} {isoforms_fpkm}")
    if transcript_gtf:
        temp_transcript_gtf = os.path.join(temp_dir, "transcripts.gtf")
        shell("mv {temp_transcript_gtf} {transcript_gtf}")
