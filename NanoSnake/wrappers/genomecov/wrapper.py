__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Get sample_id
sample = os.path.basename(snakemake.output[0]).split(".")[0]

# Run shell commands
shell("bedtools genomecov {opt} -ibam {snakemake.input} -trackopts 'type=bedGraph name={sample}' > {snakemake.output} 2> {snakemake.log}")
