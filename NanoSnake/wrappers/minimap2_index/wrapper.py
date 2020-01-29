__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
ref = snakemake.input.ref
index = snakemake.output.index

# Run shell command
shell("minimap2 -t {snakemake.threads} {opt} -d {index} {ref} &> {snakemake.log}")
