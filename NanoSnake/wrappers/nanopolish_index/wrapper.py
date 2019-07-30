__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell

# Define command depending on whether there is a summary file or not
if "summary" in snakemake.input and snakemake.input.summary:
    cmd = "nanopolish index -d {snakemake.input.fast5} {snakemake.input.fastq} -s {snakemake.input.seq_summary} -v 2> {snakemake.log}"
else:
    cmd = "nanopolish index -d {snakemake.input.fast5} {snakemake.input.fastq} -v 2> {snakemake.log}"

# Run shell command
shell (cmd)
