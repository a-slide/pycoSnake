__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.2"

# Imports
from snakemake.shell import shell

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
shell("bedtools genomecov {opt} -bg -ibam {bam} -trackopts 'type=bedGraph name={sample_id}' > {bedgraph} 2> {snakemake.log}")
