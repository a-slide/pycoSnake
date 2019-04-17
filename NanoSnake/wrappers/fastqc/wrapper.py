__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Get directory basename for fastQC
fastqc_dirbase =  os.path.dirname(snakemake.output.zip)

# Run shell command
shell("fastqc {opt} -t {snakemake.threads} --outdir {fastqc_dirbase} {snakemake.input[0]} &> {snakemake.log}")
