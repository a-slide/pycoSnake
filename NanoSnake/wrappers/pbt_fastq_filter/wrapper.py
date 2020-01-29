__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")
fastq_input = snakemake.input.fastq
fastq_output = snakemake.output.fastq

# Run shell command
shell("pyBioTools Fastq Filter -i {fastq_input} -o {fastq_output} {opt} --verbose &> {snakemake.log}")
