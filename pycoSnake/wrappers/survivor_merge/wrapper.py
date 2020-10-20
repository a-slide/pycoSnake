# Imports
from snakemake.shell import shell
import tempfile

# Wrapper info
wrapper_name = "survivor_merge"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
input_vcf = snakemake.input.vcf
output_vcf = snakemake.output.vcf

#Using temp dir for intermediate files
with tempfile.NamedTemporaryFile() as vcf_file_list:

    shell("echo '#### LS LOG ####' >> {snakemake.log}")
    shell("ls {input_vcf} > {vcf_file_list.name} 2>> {snakemake.log}")
    shell("echo '#### SURVIVOR MERGE LOG ####' >> {snakemake.log}")
    shell("SURVIVOR merge {vcf_file_list.name} {opt} {output_vcf} &>> {snakemake.log}")
