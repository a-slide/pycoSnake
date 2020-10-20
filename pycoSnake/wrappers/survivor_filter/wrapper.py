# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "survivor_filter"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
input_vcf = snakemake.input.vcf
input_bed = snakemake.input.get("bed", "NA")
output_vcf = snakemake.output.vcf

shell("SURVIVOR filter {input_vcf} {input_bed} {opt} {output_vcf} &>> {snakemake.log}")
