__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

sample_id = os.path.split(snakemake.input.bam[0])[1].rpartition(".")[0]

# Run shell command
shell ("pycoQC {opt} --verbose -f {snakemake.input.seq_summary} -a {snakemake.input.bam} --report_title {sample_id} -o {snakemake.output.html} -j {snakemake.output.json} &> {snakemake.log}")
