__author__ = "Adrien Leger"
__copyright__ = "Copyright 2019, Adrien Leger"
__email__ = "aleg@ebi.ac.uk"
__license__ = "MIT"
__version__ = "0.0.2"

# Imports
from collections import OrderedDict
import pandas as pd
from pyBioTools import Fasta
from pyfaidx import Faidx

# Shortcuts
ref_input = str(snakemake.input.ref)
ref_output = snakemake.output.ref
tsv_output = snakemake.output.get("tsv", "")

# Parse fasta file save transcript info in tabulated report and simplify transcript ids
l = []
strand_dict = {"-1":"-","1":"+"}
with open(ref_output, "w") as fa_out:

    for rec in Fasta.Reader(ref_input):
        fa_out.write(">{}\n{}\n".format(rec.short_name, rec.seq))

        # If required and possible extract all info from long sequence header
        if tsv_output:
            d = OrderedDict()
            try:
                desc = rec.long_name
                tid,_,desc = desc.partition(" ")
                d["transcript_id"]=tid
                d["length"]=len(rec)
                seqtype,_,desc = desc.partition(" ")
                d["seqtype"]=seqtype
                location,_,desc = desc.partition(" ")
                location_split = location.split(":")
                d["assembly_type"]=location_split[0]
                d["assembly_id"]=location_split[1]
                d["chrom_id"]=location_split[2]
                d["start"]=location_split[3]
                d["end"]=location_split[4]
                d["strand"]=strand_dict.get(location_split[5], ".")

                # Try to extract optional fields and store in dict
                while desc:
                    current, _,desc = desc.partition(" ")
                    k,_,v = current.partition(":")
                    if k == "description":
                        break
                    else:
                        d[k]=v
            # Skip exception silently
            except Exception:
                pass
            # Save whatever was parsed
            finally:
                if d:
                    l.append(d)

# Save tabulated report
if tsv_output:
    df = pd.DataFrame(l)
    df.to_csv(tsv_output, sep="\t", index=False)

# Index fasta file
with Faidx(ref_output) as fa_out:
    fa_out.build_index()
