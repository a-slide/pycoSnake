# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "igvtools_count"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
index = snakemake.input.ref+".fai"
bam = snakemake.input.bam
tdf = snakemake.output.tdf

# Run shell commands
shell("igvtools count {opt} {bam} {tdf} {index} &>> {snakemake.log}")
