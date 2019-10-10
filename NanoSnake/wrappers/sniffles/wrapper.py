__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
from snakemake.shell import shell
import tempfile

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Open temp directory for samtools sort temporary files
with tempfile.NamedTemporaryFile() as temp_fp:
    shell("sniffles -t {snakemake.threads} -m {snakemake.input[0]} -v {snakemake.output[0]} --tmp_file {temp_fp.name} > {snakemake.log} ")
