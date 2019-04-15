__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

from snakemake.shell import shell

if snakemake.output[0].endswith(".gz"):
    cmd = "cat {snakemake.input} | gzip -c > {snakemake.output} 2> {snakemake.log}"
else:
    cmd = "cat {snakemake.input} > {snakemake.output} 2> {snakemake.log}"

shell (cmd)
