# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "fastp"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

## Shortcuts
opt = snakemake.params.get("opt", "")
input_fastq1 = snakemake.input.fastq1
input_fastq2 = snakemake.input.get("fastq2", None)
output_fastq1 = snakemake.output.fastq1
output_fastq2 = snakemake.output.get("fastq2", None)
html = snakemake.output.html
json = snakemake.output.json

# Paired end
if input_fastq2:
    shell("fastp {opt} -w {snakemake.threads} -i {input_fastq1} -I {input_fastq2} -o {output_fastq1} -O {output_fastq2} -h {html} -j {json} &>> {snakemake.log}")

# Single end
else:
    shell("fastp {opt} -w {snakemake.threads} -i {input_fastq1} -o {output_fastq1} -h {html} -j {json} &>> {snakemake.log}")
