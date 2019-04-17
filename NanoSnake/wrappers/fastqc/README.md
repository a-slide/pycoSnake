# Wrapper for fastqc

Run Fastqc and generate an HTML report and a ZIP folder containing the results.

## Example:

```
rule fastqc:
    input:
        "reads/{sample}_reads.fastq"
    output:
        html="fastqc/{sample}_fastqc.html"
        zip="fastqc/{sample}_fastqc.zip"
    log:
        "fastqc/{sample}_fastqc.log"
    params:
        opt="--casava"
    threads:
        10
    wrapper:
        "file://path/to/wrappers/fastqc"
```
