# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "nanopolish_index"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

## Shortcuts
opt = snakemake.params.get("opt", "")
fast5 = snakemake.input.fast5
fastq = snakemake.input.fastq
seqsum = snakemake.input.get("seqsum", None)

# Run shell commands
if seqsum:
    shell("nanopolish index -v -d {fast5} {fastq} -s {seqsum} 2>> {snakemake.log}")
else:
    shell("nanopolish index -v -d {fast5} {fastq} 2>> {snakemake.log}")
