__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell
from snakemake.shell import shell

# Define command depending on whether the output has to be gzipped or not
if snakemake.output[0].endswith(".gz"):
    cmd = "cat {snakemake.input} | gzip -c > {snakemake.output} 2> {snakemake.log}"
else:
    cmd = "cat {snakemake.input} > {snakemake.output} 2> {snakemake.log}"

# Run shell command
shell (cmd)
