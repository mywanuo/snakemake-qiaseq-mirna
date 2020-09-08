# Import common module. This module defines e.g. used configuration file and 
# its validation as well as used helper functions.
include: "workflow/rules/common.smk"

# Target rule.
rule all:
    input:
        "results/counts.tsv",
        "results/counts_deduplicated.tsv"

# Pipeline modules.
include: "workflow/rules/resources.smk"
include: "workflow/rules/trimming.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/deduplication.smk"
include: "workflow/rules/counting.smk"