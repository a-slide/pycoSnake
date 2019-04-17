# Wrapper for minimap2_align

Align reads from a fastq file on a reference index. The results is then converted in bam, sorted and indexed with samtools

## Example:

```
rule minimap2_align:
    input:
        index="minimap2_index/reference.mmi"
        fastq="reads/{sample}_reads.fastq"
    output:
        "alignments/{sample}_reads.bam"
    log:
        "minimap2_align/{sample}_reads.log"
    params:
        opt="-x map-ont"
    threads:
        10
    wrapper:
        "file://path/to/wrappers/minimap2_align"
```
