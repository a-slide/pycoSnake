__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.1"

# Imports
from snakemake.shell import shell
import tempfile
import os

# Shortcuts
opt = snakemake.params.get("opt", "")
input_bam = snakemake.input.bam
input_vcf = snakemake.input.get("vcf", "") # Optional VCF to force variants
output_vcf = snakemake.output.vcf

# Run shell commands
shell("echo '#### SNIFFLES + BCFTOOLS SORT LOG ####' > {snakemake.log}")

# Using temp dir for intermediate files
with tempfile.TemporaryDirectory() as temp_dir:

    # Temorary files
    temp_vcf = os.path.join(temp_dir, "temp.vcf")
    temp_sniffles = os.path.join(temp_dir, "temp.snf")

    # Call variant with or without input vcf
    shell("echo '#### SNIFFLES LOG ####' > {snakemake.log}")
    if input_vcf:
        shell("sniffles {opt} -t {snakemake.threads} -m {input_bam} --Ivcf {input_vcf} -v {temp_vcf} --tmp_file {temp_sniffles} &>> {snakemake.log}")
    else:
        shell("sniffles {opt} -t {snakemake.threads} -m {input_bam} -v {temp_vcf} --tmp_file {temp_sniffles} &>> {snakemake.log}")

    # sort VCF
    shell("echo '#### BCFTOOLS SORT LOG ####' >> {snakemake.log}")
    shell("bcftools sort {temp_vcf} -o {output_vcf} -O v -T {temp_dir} &>> {snakemake.log}")
