__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell
import tempfile
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
threads = snakemake.threads if snakemake.threads >= 4 else 4
view_threads = sort_threads = threads//4
align_threads = threads - (view_threads + sort_threads)
outdir = os.path.dirname(os.path.abspath(snakemake.output.bam))
fastq = snakemake.input.fastq
index = snakemake.input.index
bam = snakemake.output.bam

# Run shell commands
shell("echo '#### MINIMAP2 + SAMTOOLS VIEW & SORT LOG ####' > {snakemake.log}")

with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:
    shell("minimap2 -t {align_threads} -a -L {opt} {index} {fastq} 2>> {snakemake.log}|\
        samtools view -@ {view_threads} -bh 2>> {snakemake.log} |\
        samtools sort -@ {sort_threads} -T {temp_dir} -O bam > {bam} 2>> {snakemake.log}")

shell("samtools index {bam}")
