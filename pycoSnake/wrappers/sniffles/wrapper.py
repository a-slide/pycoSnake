# Imports
from snakemake.shell import shell
import tempfile
import os

# Wrapper info
wrapper_name = "sniffles"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
input_bam = snakemake.input.bam
input_vcf = snakemake.input.get("vcf", "") # Optional VCF to force variants
output_vcf = snakemake.output.vcf

# Using temp dir for intermediate files
with tempfile.TemporaryDirectory() as temp_dir:

    # Temorary files
    temp_sniffles = os.path.join(temp_dir, "temp.snf")

    # Call variant with or without input vcf
    shell("echo '#### SNIFFLES LOG ####' >> {snakemake.log}")
    if input_vcf:
        shell("sniffles {opt} -t {snakemake.threads} -m {input_bam} --Ivcf {input_vcf} -v {output_vcf} --tmp_file {temp_sniffles} &>> {snakemake.log}")
    else:
        shell("sniffles {opt} -t {snakemake.threads} -m {input_bam} -v {output_vcf} --tmp_file {temp_sniffles} &>> {snakemake.log}")
