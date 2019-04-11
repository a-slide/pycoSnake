# An example collection of Snakemake rules imported in the main Snakefile.

# Function to get fastq path from a sample df
def get_fastq (wildcards):
    return glob (sample_df.loc[wildcards.sample, "fastq"])
