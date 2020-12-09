# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "nanopolish_call_methylation"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
fastq = snakemake.input.fastq
bam = snakemake.input.bam
ref = snakemake.input.ref
tsv = snakemake.output.tsv

# Run shell commands
shell("nanopolish call-methylation {opt} -t {snakemake.threads} -r {fastq} -b {bam} -g {ref} > {tsv} 2>> {snakemake.log}")
