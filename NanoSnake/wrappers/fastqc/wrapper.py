__author__ = "Adrien Leger & Julian de Ruiter"
__copyright__ = "Copyright 2019, Adrien Leger & Julian de Ruiter"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"

from os import path
from tempfile import TemporaryDirectory
import shutil
from snakemake.shell import shell

def basename_without_ext(fn):
    """Returns basename of file path, without the file extension."""
    split_ind = 2 if fn.endswith(".gz") else 1
    base = ".".join(path.basename(fn).split(".")[:-split_ind])
    return base

# Get optional args if unavailable
opt = snakemake.params.get("opt", "")

# Run fastqc, since there can be race conditions if multiple jobs use the same fastqc dir, we create a temp dir.
with TemporaryDirectory() as tempdir:
    shell("fastqc {opt} -t {snakemake.threads} --outdir {tempdir} {snakemake.input[0]} &> {snakemake.log}")

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir, output_base + "_fastqc.html")
    _ = shutil.move(html_path, snakemake.output.html)
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")
    _ = shutil.move(zip_path, snakemake.output.zip)
