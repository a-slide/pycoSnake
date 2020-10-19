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
input_vcf = snakemake.input.vcf
output_vcf = snakemake.output.vcf

#Using temp dir for intermediate files
with tempfile.NamedTemporaryFile() as vcf_file_list:

    shell("echo '#### LS LOG ####' > {snakemake.log}")
    shell("ls {input_vcf} > {vcf_file_list.name} 2>> {snakemake.log}")
    shell("cat {vcf_file_list.name}")
    shell("echo '#### SURVIVOR MERGE LOG ####' >> {snakemake.log}")
    shell("SURVIVOR merge {vcf_file_list.name} {opt} {output_vcf} &>> {snakemake.log}")
