__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run shell command
shell ("pycoQC {opt} --verbose -f {snakemake.input.seq_summary} -a {snakemake.input.bam} -o {snakemake.output.html} -j {snakemake.output.json} &> {snakemake.log}")
