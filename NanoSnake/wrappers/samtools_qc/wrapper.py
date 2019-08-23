__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Run shell command
shell("echo '#### SAMTOOLS STATS LOG ####' > {snakemake.log}")
shell("samtools stats {snakemake.input[0]} > {snakemake.output.stats} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS FLAGSTAT LOG ####' >> {snakemake.log}")
shell("samtools flagstat {snakemake.input[0]} > {snakemake.output.flagstat} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS IDXSTATS LOG ####' >> {snakemake.log}")
shell("samtools idxstats {snakemake.input[0]} > {snakemake.output.idxstats} 2>> {snakemake.log}")
