__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.2"

# Imports
from pyBioTools import Fasta
from pyfaidx import Faidx

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
