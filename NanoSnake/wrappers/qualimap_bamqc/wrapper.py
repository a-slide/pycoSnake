__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

# Imports
import os
from snakemake.shell import shell

shell("unset DISPLAY && qualimap bamqc --java-mem-size={snakemake.resources.mem_mb}M -bam {snakemake.input[0]} -outdir {snakemake.output[0]} > {snakemake.log}")
