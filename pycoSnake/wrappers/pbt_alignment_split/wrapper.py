# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "pbt_alignment_split"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
bam_input = snakemake.input.bam
bam_output = snakemake.output.bam

# Run shell command
shell("pyBioTools --version >> {snakemake.log}")
shell("pyBioTools Alignment Split {opt} -i {bam_input} -l {bam_output} --verbose &>> {snakemake.log}")
