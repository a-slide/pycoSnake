__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

import textwrap
from snakemake.shell import shell

# Get optional args if unavailable
if "stats" in snakemake.output:
    shell("samtools stats {snakemake.input[0]} > ")
if "flagstat" in snakemake.output:
    shell("samtools flagstat {snakemake.input[0]} > ")
if "idxstats" in snakemake.output:
    shell("samtools idxstats {snakemake.input[0]} > ")
