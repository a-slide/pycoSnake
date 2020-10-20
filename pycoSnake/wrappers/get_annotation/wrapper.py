# Imports
from snakemake.shell import shell
import tempfile
import os

# Wrapper info
wrapper_name = "get_annotation"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
input_gff3 = str(snakemake.input.gff3)
output_gff3 = snakemake.output.gff3
output_gtf = snakemake.output.gtf
outdir = os.path.dirname(os.path.abspath(output_gff3))

with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:

    # Gzipped files
    if input_gff3.rpartition(".")[-1].lower()=="gz":
        temp_gff3 = os.path.join(temp_dir, "temp.gff3")
        shell("gunzip -c {input_gff3} > {temp_gff3}")
        input_gff3 = temp_gff3

    # Convert and clean annotations
    shell("gffread {input_gff3} {opt} -F --keep-genes -o {output_gff3} &>> {snakemake.log}")
    shell("gffread {input_gff3} {opt} -F --keep-genes -T -o {output_gtf} &>> {snakemake.log}")
