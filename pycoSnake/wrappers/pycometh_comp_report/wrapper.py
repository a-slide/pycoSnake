# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "pycometh_comp_report"
wrapper_version = "0.0.4"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
tsv = snakemake.input.tsv
gff3 = snakemake.input.gff3
ref = snakemake.input.ref
outdir = snakemake.output.outdir

# Run shell command
shell("pycoMeth --version >> {snakemake.log}")
shell(f"pycoMeth Comp_Report {opt} -i {tsv} -f {ref} -g {gff3} -o {outdir} 2>> {snakemake.log}")
