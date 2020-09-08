import pandas as pd
from snakemake.utils import min_version
from snakemake.utils import validate

# Define minimum Snakemake version.
min_version("5.14.0")

# Define and validate configuration file.
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Import and validate sample list.
samples = pd.read_table("samples.tsv")
validate(samples, schema="../schemas/samples.schema.yaml")

# Define helper functions that extract wildcards from sample list and return 
# them when required.

def get_wildcards():
    """Extracts wildcards from sample list and returns them as list."""
    return samples["sample"].unique().tolist()

def get_fastq(wildcards):
    """Returns list of .fastq files located at the /fastq directory."""

    # Extract sample specific .fastq files from sample sheet.
    fastqs = samples[samples["sample"] == wildcards.sample]["fastq"].tolist()
    
    # Return list of .fastq files.
    return fastqs