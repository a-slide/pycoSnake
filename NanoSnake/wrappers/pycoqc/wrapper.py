__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
import os
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
seqsum = snakemake.input.seqsum
bam = snakemake.input.get("bam", None)
html = snakemake.output.get("html", None)
json = snakemake.output.get("json", None)

# Get sample_id
sample_id = snakemake.params.get("sample_id", None)
if not sample_id:
    try:
        sample_id = os.path.split(bam)[1].rpartition(".")[0]
    except Exception:
        sample_id = "Sample"

# Define conditional inputs
if bam:
    input = f"-f {seqsum} -a {bam}"
else:
    input = f"-f {seqsum}"

# Define conditional outputs
if html and json:
    output = f"--report_title {sample_id} -o {html} -j {json}"
elif html:
    output = f"--report_title {sample_id} -o {html}"
elif json:
    output = f"-j {json}"

# Run shell command
shell (f"pycoQC {opt} {input} {output} &> {snakemake.log}")
