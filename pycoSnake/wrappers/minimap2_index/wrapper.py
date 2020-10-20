# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "minimap2_index"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
ref = snakemake.input.ref
index = snakemake.output.index

# Run shell command
shell("minimap2 -t {snakemake.threads} {opt} -d {index} {ref} &>> {snakemake.log}")
