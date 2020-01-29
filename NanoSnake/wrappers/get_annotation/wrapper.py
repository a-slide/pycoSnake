__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.2"

# Imports
from snakemake.shell import shell
import tempfile
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
input_gff3 = str(snakemake.input.gff3)
output_gff3 = snakemake.output.gff3
output_gtf = snakemake.output.gtf
outdir = os.path.dirname(os.path.abspath(output_gff3))

# Create temp dir
with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:
    # Gzipped files
    if input_gff3.rpartition(".")[-1].lower()=="gz":
        temp_gff3 = os.path.join(temp_dir, "temp.gff3")
        shell("gunzip -c {input_gff3} > {temp_gff3}")
        input_gff3 = temp_gff3

    # Convert and clean annotations
    shell("gffread {input_gff3} {opt} -F -o {output_gff3} &>> {snakemake.log}")
    shell("gffread {input_gff3} {opt} -F -T -o {output_gtf} &>> {snakemake.log}")
