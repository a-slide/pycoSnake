__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell
import tempfile
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
bam_input = snakemake.input.bam
bam_output = snakemake.output.bam

# Run shell command
shell("pyBioTools Alignment Filter -i {bam_input} -o {bam_output} {opt} --verbose &> {snakemake.log}")