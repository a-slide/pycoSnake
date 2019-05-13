__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell
import tempfile

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Threads sharing
threads = snakemake.threads if snakemake.threads >= 4 else 4
view_threads = sort_threads = threads//4
align_threads = threads - (view_threads + sort_threads)

# Run shell commands
shell("echo '#### MINIMAP2 + SAMTOOLS VIEW & SORT LOG ####' > {snakemake.log}")

# Open temp directory for samtools sort temporary files
with tempfile.TemporaryDirectory() as temp_dir:
    shell("minimap2 -t {align_threads} -a -L {opt} {snakemake.input.index} {snakemake.input.fastq} 2>> {snakemake.log}|\
        samtools view -@ {view_threads}  -bh 2>> {snakemake.log} |\
        samtools sort -@ {sort_threads} -T {temp_dir} -O bam > {snakemake.output[0]} 2>> {snakemake.log}")

shell("samtools index {snakemake.output[0]}")
