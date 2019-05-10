# Wrapper for bedtools genomecov



## Example:

```
rule genomecov:
    input:
        rules.merge_fastq.output
    output:
        path.join("results", genomecov_dir,"{sample}.bedgraph")
    log:
        path.join("logs", genomecov_dir,"{sample}.log")
    params:
        opt=config["genomecov"]["opt"],
    wrapper:
        "genomecov"
```
