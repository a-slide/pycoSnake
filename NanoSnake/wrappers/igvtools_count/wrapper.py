__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run shell commands
shell("touch {snakemake.log}")

# Create ref genome index if needed
index = snakemake.input.ref+".fai"
if not os.access(index, os.R_OK):
    shell("echo '#### SAMTOOLS INDEX LOG####' >> {snakemake.log}")
    shell("samtools faidx {snakemake.input.ref} 2>> {snakemake.log}")

shell("echo '#### IGVTOOLS COUNT LOG ####' &>> {snakemake.log}")
shell("igvtools count {opt} {snakemake.input.bam} {snakemake.output} {index} &>> {snakemake.log}")
