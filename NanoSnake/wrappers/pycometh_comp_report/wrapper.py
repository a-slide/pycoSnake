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
tsv = snakemake.input.tsv
gff3 = snakemake.input.gff3
outdir = snakemake.output.outdir

# Run shell command
shell(f"pycoMeth Comp_Report {opt} -i {tsv} -g {gff3} -o {outdir} 2> {snakemake.log}")
