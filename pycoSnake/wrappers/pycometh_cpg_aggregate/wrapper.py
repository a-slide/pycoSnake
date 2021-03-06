# Imports
from snakemake.shell import shell
import os

# Wrapper info
wrapper_name = "pycometh_cpg_aggregate"
wrapper_version = "0.0.5"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
input_tsv = snakemake.input.tsv
ref = snakemake.input.ref
output_tsv = snakemake.output.get("tsv", "")
output_bed = snakemake.output.get("bed", "")
output_bed_index = snakemake.output.get("bed_index", "")

# Define conditional IO
output = ""
if output_bed:
    sample_id = snakemake.params.get("sample_id", None)
    if not sample_id:
        sample_id = os.path.basename(output_bed).rpartition(".")[0]
    output += f" -b {output_bed} -s {sample_id} "
if output_tsv:
    output += f" -t {output_tsv} "

# Run shell command
shell("pycoMeth --version >> {snakemake.log}")
shell(f"pycoMeth CpG_Aggregate {opt} -i {input_tsv} -f {ref} {output} 2>> {snakemake.log}")

# Optional Indexing with igv
if output_bed_index and output_bed_index.endswith(".idx") and output_bed and output_bed.endswith(".bed"):
    shell(f"igvtools index {output_bed} &>> {snakemake.log}")
