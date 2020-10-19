__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.3"

# Imports
import os
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
input_tsv = snakemake.input.tsv
ref = snakemake.input.ref
output_tsv = snakemake.output.get("tsv", "")
output_bed = snakemake.output.get("bed", "")
output_bed_index = snakemake.output.get("bed_index", "")

# Get sample_id
sample_id = snakemake.params.get("sample_id", None)
if not sample_id:
    try:
        if input_tsv.endswith(".gz"):
            sample_id = os.path.basename(input_tsv).rpartition(".")[0].rpartition(".")[0]
        else:
            sample_id = os.path.basename(input_tsv).rpartition(".")[0]
    except Exception:
        sample_id = "Sample"

# Define conditional IO
output = ""
if output_tsv:
    output += f" -t {output_tsv} "
if output_bed:
    output += f" -b {output_bed}"

# Run shell command
shell(f"pycoMeth CpG_Aggregate {opt} -i {input_tsv} -f {ref} {output} -s {sample_id} 2> {snakemake.log}")

# Optional Indexing with igv
if output_bed_index and output_bed_index.endswith(".idx") and output_bed and output_bed.endswith(".bed"):
    shell(f"igvtools index {output_bed} &>> {snakemake.log}")
