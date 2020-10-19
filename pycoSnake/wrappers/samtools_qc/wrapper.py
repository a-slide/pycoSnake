__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Shortcuts
bam = snakemake.input.bam
stats = snakemake.output.stats
flagstat = snakemake.output.flagstat
idxstats = snakemake.output.idxstats

# Run shell command
shell("echo '#### SAMTOOLS STATS LOG ####' > {snakemake.log}")
shell("samtools stats {bam} > {stats} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS FLAGSTAT LOG ####' >> {snakemake.log}")
shell("samtools flagstat {bam} > {flagstat} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS IDXSTATS LOG ####' >> {snakemake.log}")
shell("samtools idxstats {bam} > {idxstats} 2>> {snakemake.log}")
