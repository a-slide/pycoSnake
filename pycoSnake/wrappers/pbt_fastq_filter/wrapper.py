# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "pbt_fastq_filter"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")
fastq_input = snakemake.input.fastq
fastq_output = snakemake.output.fastq

# Run shell command
shell("pyBioTools --version >> {snakemake.log}")
shell("pyBioTools Fastq Filter -i {fastq_input} -o {fastq_output} {opt} --verbose &>> {snakemake.log}")
