# Imports
from snakemake.shell import shell
import os

# Wrapper info
wrapper_name = "pycometh_comp_report"
wrapper_version = "0.0.5"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
tsv = snakemake.input.tsv
gff3 = snakemake.input.gff3
ref = snakemake.input.ref
summary_report = snakemake.output.summary_report
outdir=os.path.dirname(summary_report)

# Run shell command
shell("pycoMeth --version >> {snakemake.log}")
shell(f"pycoMeth Comp_Report {opt} -i {tsv} -f {ref} -g {gff3} -o {outdir} 2>> {snakemake.log}")
