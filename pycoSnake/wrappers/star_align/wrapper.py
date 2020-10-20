# Imports
from snakemake.shell import shell
import tempfile
import os

# Wrapper info
wrapper_name = "star_align"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
opt = snakemake.params.get("opt", "")
index_dir = os.path.abspath(snakemake.input.index_dir)+"/"
fastq1 = snakemake.input.fastq1
fastq2 = snakemake.input.fastq2
sj = snakemake.output.get("sj", None)
count = snakemake.output.get("count", None)
bam = snakemake.output.get("bam", None)
bam_index = snakemake.output.get("bam_index", None)
star_log = snakemake.output.get("star_log", None)

# Unzipping option
if fastq1.endswith(".gz") and fastq2.endswith(".gz"):
    unzip_option = "--readFilesCommand 'gunzip -c'"
elif fastq1.endswith(".gz") or fastq2.endswith(".gz"):
    raise ValueError ("both fastq files needs to be either gziped or uncompressed")
else:
    unzip_option = ""

# Run shell command
outdir = os.path.dirname(os.path.abspath(snakemake.output[0]))
with tempfile.TemporaryDirectory(dir=outdir) as temp_dir:
    temp_dir = temp_dir+os.path.sep

    # Run shell command
    shell("STAR {opt} {unzip_option}\
        --runMode alignReads\
        --outSAMtype BAM SortedByCoordinate\
        --runThreadN {snakemake.threads}\
        --genomeDir {index_dir}\
        --readFilesIn {fastq1} {fastq2}\
        --outFileNamePrefix {temp_dir}\
        --quantMode GeneCounts\
        &>> {snakemake.log}")

    if sj:
        temp_file = os.path.join(temp_dir, "SJ.out.tab")
        shell("mv {temp_file} {sj}")
    if count:
        temp_file = os.path.join(temp_dir, "ReadsPerGene.out.tab")
        shell("mv {temp_file} {count}")
    if bam:
        temp_file = os.path.join(temp_dir, "Aligned.sortedByCoord.out.bam")
        shell("mv {temp_file} {bam}")
        if bam_index:
            shell("samtools index {bam}")
    if star_log:
        temp_file = os.path.join(temp_dir, "Log.final.out")
        shell("mv {temp_file} {star_log}")
