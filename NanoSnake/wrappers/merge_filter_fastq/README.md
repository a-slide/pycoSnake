# Wrapper for merge_filter_fastq.

Concatenate multiple fastq into a single one and optionnally filters reads based on quality and seq len
If the output has a ".gz extension" the output is automatically compressed in gzip format

## Example:

```
rule merge_filter_fastq:
    input:
        "dir_containing_fastq"
    output:
        "reads/{sample}_merged.fastq"
    log:
        "logs/concat_files/{sample}.log"
    wrapper:
        "file://path/to/wrappers/merge_filter_fastq"
```
