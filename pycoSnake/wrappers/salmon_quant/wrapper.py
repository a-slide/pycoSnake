# Imports
from snakemake.shell import shell
import os

# Wrapper info
wrapper_name = "salmon_quant"
wrapper_version = "0.0.4"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
index_dir = snakemake.input.index_dir
fastq1 = snakemake.input.fastq1
fastq2 = snakemake.input.fastq2
quant_dir = snakemake.output.quant_dir
counts = snakemake.output.get("counts", None)
os.makedirs(quant_dir, exist_ok=True)

# Run shell command
shell("salmon quant {opt} -p {snakemake.threads} -i {index_dir} -1 {fastq1} -2 {fastq2} -o {quant_dir} &>> {snakemake.log}")

# Move counts if given
if counts:
    src = os.path.join(quant_dir, "quant.sf")
    shell("mv {src} {counts}")
