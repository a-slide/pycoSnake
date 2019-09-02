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
shell("touch {snakemake.log}")

# Run shell commands
shell("echo '#### NANOPOLISH INDEX LOG ####' >> {snakemake.log}")
shell("nanopolish index -d {snakemake.input.fast5} {snakemake.input.fastq} -s {snakemake.input.seq_summary} -v 2>> {snakemake.log}")

shell("echo '#### NANOPOLISH CALL-METHYLATION LOG ####' > {snakemake.log}")
shell("nanopolish call-methylation -t {snakemake.threads} {opt} -r {snakemake.input.fastq} -b {snakemake.input.bam} -g {snakemake.input.ref} > {snakemake.output} 2>> {snakemake.log}")
