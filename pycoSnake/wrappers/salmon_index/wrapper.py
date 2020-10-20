# Imports
from snakemake.shell import shell
import os

# Wrapper info
wrapper_name = "salmon_index"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
ref = snakemake.input.ref
index_dir = snakemake.output.index_dir
os.makedirs(index_dir, exist_ok=True)

# Run shell command
shell("salmon index {opt} -p {snakemake.threads} -t {ref} -i {index_dir} &>> {snakemake.log}")
