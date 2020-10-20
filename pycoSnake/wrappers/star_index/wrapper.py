# Imports
from snakemake.shell import shell
from pyfaidx import Fasta
from math import log2
import os

# Wrapper info
wrapper_name = "star_index"
wrapper_version = "0.0.4"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
ref = snakemake.input.ref
annotation = snakemake.input.annotation
index_dir = os.path.abspath(snakemake.output.index_dir)+"/"
os.makedirs(index_dir, exist_ok=True)

# Comput index base depending on genome length
genome_len = 0
with Fasta(ref) as fa:
    for seq in fa:
        genome_len+=len(seq)
indexNbases = min(14, int(log2(genome_len)/2) - 1)

# Run shell command
shell("STAR {opt} \
    --genomeSAindexNbases {indexNbases} \
    --runMode genomeGenerate \
    --runThreadN {snakemake.threads} \
    --genomeDir {index_dir} \
    --genomeFastaFiles {ref} \
    --sjdbGTFfile {annotation} \
    --outFileNamePrefix {index_dir} \
    &>> {snakemake.log}")
