# -*- coding: utf-8 -*-

# Copyright 2019 Adrien Leger
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd
from snakemake.utils import validate, min_version

# set minimum snakemake version
min_version("5.4.2")

#~~~~~~~~~~~~~~~~~~~~~~~~load config and sample sheets~~~~~~~~~~~~~~~~~~~~~~~~#

if config:
    print (config)
else:
    configfile: "config.yaml"

sample_df = pd.read_csv (config["sample_fn"], sep="\t", index_col=0)

wildcard_constraints:
    sample="|".join(sample_df.index),
    fastq=list(sample_df["fastq"]),
    condition="|".join(sample_df["condition"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Rules~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule samples:
    input: exp
    output: "test_output.fa"
    shell: "cp {input} {output}"

rule r2:
    input: "test_output.fa"
    shell: "cat {input}"
