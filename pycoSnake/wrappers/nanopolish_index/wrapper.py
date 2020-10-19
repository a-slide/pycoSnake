__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
import os
from snakemake.shell import shell

## Shortcuts
opt = snakemake.params.get("opt", "")
fast5 = snakemake.input.fast5
fastq = snakemake.input.fastq
seqsum = snakemake.input.get("seqsum", None)

# Run shell commands
if seqsum:
    shell("nanopolish index -v -d {fast5} {fastq} -s {seqsum} 2> {snakemake.log}")
else:
    shell("nanopolish index -v -d {fast5} {fastq} 2> {snakemake.log}")
