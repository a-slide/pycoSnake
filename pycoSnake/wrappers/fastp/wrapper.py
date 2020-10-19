__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

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
    shell("fastp {opt} -w {snakemake.threads} -i {input_fastq1} -I {input_fastq2} -o {output_fastq1} -O {output_fastq2} -h {html} -j {json} &> {snakemake.log}")

# Single end
else:
    shell("fastp {opt} -w {snakemake.threads} -i {input_fastq1} -o {output_fastq1} -h {html} -j {json} &> {snakemake.log}")
