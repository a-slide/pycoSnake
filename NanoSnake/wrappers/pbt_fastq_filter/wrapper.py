__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run shell command
shell("pyBioTools Fastq Filter -i {snakemake.input[0]} -o {snakemake.output[0]} {opt} --verbose &> {snakemake.log}")
