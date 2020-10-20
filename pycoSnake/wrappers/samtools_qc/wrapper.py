# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "samtools_qc"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
bam = snakemake.input.bam
stats = snakemake.output.stats
flagstat = snakemake.output.flagstat
idxstats = snakemake.output.idxstats

# Run shell command
shell("echo '#### SAMTOOLS STATS LOG ####' >> {snakemake.log}")
shell("samtools stats {bam} > {stats} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS FLAGSTAT LOG ####' >> {snakemake.log}")
shell("samtools flagstat {bam} > {flagstat} 2>> {snakemake.log}")

shell("echo '#### SAMTOOLS IDXSTATS LOG ####' >> {snakemake.log}")
shell("samtools idxstats {bam} > {idxstats} 2>> {snakemake.log}")
