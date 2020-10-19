__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
bam = snakemake.input.bam
gtf = snakemake.input.gtf
counts = snakemake.output.counts

# Run shell command
shell("featureCounts {opt} -T {snakemake.threads} -a {gtf} -o {counts} {bam} &> {snakemake.log}")
