# Wrapper for minimap2_index

Generate an index for minimap2 from a fasta reference genome

## Example:

```
rule minimap2_index:
    input:
        "reference/reference.fa"
    output:
        "minimap2_index/reference.mmi"
    log:
        "minimap2_index/reference.log"
    params:
        ""
    threads:
        10
    wrapper:
        "file://path/to/wrappers/minimap2/index"
```
