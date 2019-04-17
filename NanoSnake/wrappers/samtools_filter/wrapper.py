__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Split threads in 2 for view ans sort
threads = 1 if snakemake.threads <= 1 else snakemake.threads//2

# Run shell command
shell("samtools view -bh {opt} -@ {threads} {snakemake.input[0]} | samtools sort -O bam > {snakemake.output[0]}  2> {snakemake.log}")
shell("samtools index {snakemake.output[0]}")
