__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
import os
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
fastq = snakemake.input.fastq
bam = snakemake.input.bam
ref = snakemake.input.ref
tsv = snakemake.output.tsv

# Run shell commands
shell("nanopolish call-methylation {opt} -t {snakemake.threads} -r {fastq} -b {bam} -g {ref} > {tsv} 2> {snakemake.log}")
