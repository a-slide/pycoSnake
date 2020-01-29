__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell
import tempfile

# Shortcuts
opt = snakemake.params.get("opt", "")
bam = snakemake.input.bam
vcf = snakemake.output.vcf

# Open temp directory for samtools sort temporary files
with tempfile.NamedTemporaryFile() as temp_fp:
    shell("sniffles -t {snakemake.threads} -m {bam} -v {vcf} --tmp_file {temp_fp.name} > {snakemake.log}")
