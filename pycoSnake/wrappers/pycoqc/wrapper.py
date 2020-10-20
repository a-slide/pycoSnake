# Imports
from snakemake.shell import shell
import os

# Wrapper info
wrapper_name = "pycoqc"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

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
shell(f"pycoQC --version >> {snakemake.log}")
shell (f"pycoQC {opt} {input} {output} &>> {snakemake.log}")
