# Imports
import gzip

# Wrapper info
wrapper_name = "nanopolish_concat"
wrapper_version = "0.0.1"
author = "Adrien Leger"
license = "MIT"

# Shortcuts
input_tsv_list = snakemake.input.tsv_list
output_tsv = snakemake.output.tsv

with open (str(snakemake.log), "w") as log_fp:
    log_fp.write(f'Wrapper {wrapper_name} v{wrapper_version} / {author} / Licence {license}\n')

    # Concatenate and manage gzip compression
    open_fun, open_mode = (gzip.open, "wt") if output_tsv.endswith(".gz") else (open, "w")

    with open_fun(output_tsv, open_mode) as output_fp:
        first = True
        for input_tsv in input_tsv_list:
            log_fp.write(f'Reading file {input_tsv}\n')
            open_fun, open_mode = (gzip.open, "rt") if input_tsv.endswith(".gz") else (open, "r")
            with open_fun(input_tsv, open_mode) as input_fp:
                if not first:
                    # flush header if not first file
                    _ = input_fp.readline()
                else:
                    first=False
                output_fp.write(input_fp.read())
