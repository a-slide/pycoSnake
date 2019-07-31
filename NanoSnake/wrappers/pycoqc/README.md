# Wrapper for fastqc

Run pycoQC and generate an interactive HTML and a text json report

## Example:

```
rule fastqc:
    input:
        summary="data/{sample}/sequencing_summary.txt"
    output:
        html="pycoQC/{sample}.html",
        json="pycoQC/{sample}.json"
    log:
        "pycoQC/{sample}.log"
    params:
        opt=""
    wrapper:
        "file://path/to/wrappers/pycoqc"
```
