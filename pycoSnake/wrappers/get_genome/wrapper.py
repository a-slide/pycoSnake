# Imports
from snakemake.shell import shell
from pyBioTools import Fasta
from pyfaidx import Faidx

# Wrapper info
wrapper_name = "get_genome"
wrapper_version = "0.0.3"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
ref_input = str(snakemake.input.ref)
ref_output = snakemake.output.ref

# Parse fasta file save transcript info in tabulated report and simplify transcript ids
with open(ref_output, "w") as fa_out:
    for rec in Fasta.Reader(ref_input):
        fa_out.write(">{}\n{}\n".format(rec.short_name, rec.seq))

# Index fasta file
with Faidx(ref_output) as fa_out:
    fa_out.build_index()
