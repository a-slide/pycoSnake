# Wrapper for for multiple QC on aligned bam files

This wrapper runs a series of QC on the bam files including qualimap_bamqc, samtools stats, samtools flagstat and samtools idxstats.
The following layout works well with MultiQC

## Example:

```
rule bamqc:
    input:
        reads_{sample}.bam
    output:
        qualimap="bamqc/{sample}/qualimapReport.html",
        stats="bamqc/{sample}_samtools_stats.txt",
        flagstat="bamqc/{sample}_samtools_flagstat.txt",
        idxstats="bamqc/{sample}_samtools_idxstats.txt",
    log:
        "bamqc/{sample}.log,
    wrapper:
        "file://path/to/wrappers/bamqc"
```
