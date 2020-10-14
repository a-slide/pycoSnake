__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell

# Shortcuts
opt = snakemake.params.get("opt", "")
input_vcf = snakemake.input.vcf
input_bed = snakemake.input.get("bed", "NA")
output_vcf = snakemake.output.vcf

shell("SURVIVOR filter {input_vcf} {input_bed} {opt} {output_vcf} &> {snakemake.log}")
