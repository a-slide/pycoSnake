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
ref = snakemake.input.ref
output_tsv = snakemake.output.get("tsv", "")
output_bed = snakemake.output.get("bed", "")
output_bed_index = snakemake.output.get("bed_index", "")

# Define conditional IO
output = ""
if output_tsv:
    output += f" -t {output_tsv} "
if output_bed:
    output += f" -b {output_bed}"

# Run shell command
shell(f"pycoMeth CGI_Finder {opt} -f {ref} {output} 2> {snakemake.log}")

# Optional Indexing with igv
if output_bed_index and output_bed_index.endswith(".idx") and output_bed and output_bed.endswith(".bed"):
    shell(f"igvtools index {output_bed} &>> {snakemake.log}")
