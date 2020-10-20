# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "subread_featurecounts"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Imports
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
bam = snakemake.input.bam
gtf = snakemake.input.gtf
counts = snakemake.output.counts

# Run shell command
shell("featureCounts {opt} -T {snakemake.threads} -a {gtf} -o {counts} {bam} &>> {snakemake.log}")
