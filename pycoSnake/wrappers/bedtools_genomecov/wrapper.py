# Imports
from snakemake.shell import shell

# Wrapper info
wrapper_name = "bedtools_genomecov"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
bam = snakemake.input.bam
bedgraph = snakemake.output.bedgraph

# Get sample_id
sample_id = snakemake.params.get("sample_id", None)
if not sample_id:
    try:
        sample_id = os.path.split(bam)[1].rpartition(".")[0]
    except Exception:
        sample_id = "Sample"

# Run shell commands
shell("bedtools genomecov {opt} -bg -ibam {bam} -trackopts 'type=bedGraph name={sample_id}' > {bedgraph} 2>> {snakemake.log}")
