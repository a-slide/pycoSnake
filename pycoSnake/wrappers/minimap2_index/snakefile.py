# Imports
from os.path import join

# Input and output data
ref = join(config["data_dir"], "reference", "ref.fa")
index = "index/ref.mmi"

# Rules
rule all:
    input: [index]

rule minimap2_index:
    input: ref=ref
    output: index=index
    threads: 2
    params: opt=""
    resources: mem_mb=1000
    log: "minimap2_index.log"
    wrapper: "minimap2_index"
