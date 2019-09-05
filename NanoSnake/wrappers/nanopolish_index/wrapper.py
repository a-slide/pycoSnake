__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run shell commands
shell("nanopolish index -d {snakemake.input.fast5} {snakemake.input.fastq} -s {snakemake.input.seq_summary} -v 2> {snakemake.log}")
