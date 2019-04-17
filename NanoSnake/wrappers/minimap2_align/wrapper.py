__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run shell commands
shell("""
    minimap2 -t {snakemake.threads} -a -L {opt} {snakemake.input.index} {snakemake.input.fastq} 2> {snakemake.log} |
    samtools view -bh | samtools sort -O bam > {snakemake.output[0]}
    """)
shell("samtools index {snakemake.output[0]}")
