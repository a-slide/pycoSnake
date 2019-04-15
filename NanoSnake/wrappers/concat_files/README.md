# Wrapper for concat_files.

Concatenate multiple file into a single one.
If the output has a ".gz extension" the output is automatically compressed in gzip format

## Example:

```
rule concat_files:
    input:
        ["reads/{sample}_1.fastq", "reads/{sample}_2.fastq"]
    output:
        "reads/{sample}_merged.fastq"
    log:
        "logs/concat_files/{sample}.log"
    wrapper:
        "file://path/to/wrappers/concat_files"
```
