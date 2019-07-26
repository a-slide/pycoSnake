__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell

# Create empty files
shell("touch {snakemake.output}")
shell("touch {snakemake.log}")

# For some reason if the list of files is too long it raises an error if using cat directly
for file in snakemake.input:
    # Define command depending on whether the output has to be gzipped or not
    if snakemake.output[0].endswith(".gz"):
        shell("gzip -c {file} >> {snakemake.output} 2>> {snakemake.log}")
    else:
        shell("cat {file} >> {snakemake.output} 2>> {snakemake.log}")
