__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

from snakemake.shell import shell

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")
samtools_view_opt = snakemake.params.get("samtools_view_opt", "")

shell("""
    minimap2 -t {snakemake.threads} -a -L {opt} {snakemake.input.index} {snakemake.input.fastq} |
    samtools view -bh | samtools sort -O bam > {snakemake.output[0]} 2> {snakemake.log}
    """)

shell("samtools index {snakemake.output[0]}")
