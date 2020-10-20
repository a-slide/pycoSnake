# Imports
from snakemake.shell import shell
import tempfile
import os

# Wrapper info
wrapper_name = "ngmlr"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
threads = snakemake.threads if snakemake.threads >= 4 else 4
view_threads = sort_threads = threads//4
align_threads = threads - (view_threads + sort_threads)
outdir = os.path.dirname(os.path.abspath(snakemake.output.bam))
fastq = snakemake.input.fastq
ref = snakemake.input.ref
bam = snakemake.output.bam

# Run shell commands
shell("echo '#### NGMLR + SAMTOOLS LOG ####' >> {snakemake.log}")

with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:
    shell("ngmlr -t {align_threads} {opt} -r {ref} -q {fastq} 2>> {snakemake.log}|\
        samtools view -@ {view_threads}  -bh 2>> {snakemake.log} |\
        samtools sort -@ {sort_threads} -T {temp_dir} -O bam > {bam} 2>> {snakemake.log}")

shell("samtools index {bam}")
