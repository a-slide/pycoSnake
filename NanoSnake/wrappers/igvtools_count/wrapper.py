__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
index = snakemake.input.ref+".fai"
bam = snakemake.input.bam
tdf = snakemake.output.tdf

# Run shell commands
shell("igvtools count {opt} {bam} {tdf} {index} &> {snakemake.log}")
