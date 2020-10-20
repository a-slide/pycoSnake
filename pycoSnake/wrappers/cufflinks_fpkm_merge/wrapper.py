# Imports
from snakemake.shell import shell
import pandas as pd
import os

# Wrapper info
wrapper_name = "cufflinks_fpkm_merge"
wrapper_version = "0.0.2"
author = "Adrien Leger"
license = "MIT"
shell("echo 'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}' > {snakemake.log}")

# Shortcuts
input_fpkm_genes = snakemake.input.get("fpkm_genes", None)
input_fpkm_isoforms = snakemake.input.get("fpkm_isoforms", None)
output_fpkm_genes = snakemake.output.get("fpkm_genes", None)
output_fpkm_isoforms = snakemake.output.get("fpkm_isoforms", None)

for input_fpkm, output_fpkm in ((input_fpkm_genes,output_fpkm_genes),(input_fpkm_isoforms,output_fpkm_isoforms)):
    if input_fpkm and output_fpkm:
        df = pd.DataFrame()
        for c in input_fpkm:
            sample_id = os.path.basename(c).rpartition(".")[0]
            sample_df = pd.read_csv(c, sep="\t")
            sample_df = sample_df[["tracking_id", "FPKM"]]
            sample_df = sample_df.dropna()
            sample_df = sample_df.rename(columns={"tracking_id":"id","FPKM":sample_id})
            sample_df = sample_df.set_index("id")
            df = df.join(sample_df, how='outer')

        df = df.fillna(0)
        df.to_csv(output_fpkm, sep="\t")
