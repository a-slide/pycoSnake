__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.3"

# Imports
from snakemake.shell import shell
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
ref = snakemake.input.ref
index_dir = snakemake.output.index_dir
os.makedirs(index_dir, exist_ok=True)

# Run shell command
shell("salmon index {opt} -p {snakemake.threads} -t {ref} -i {index_dir} &> {snakemake.log}")
