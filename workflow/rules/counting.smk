# Count raw read counts from BAM file and extract relevant columns from output.
# Results will be written to tab-delimited .tsv file.
rule count_raw_read_counts:
    input:
        bam="results/{sample}/bam/{sample}_merged.bam",
        bai="results/{sample}/bam/{sample}_merged.bai"
    output:
        counts="results/{sample}/{sample}_counts.tsv"
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/count_raw_read_counts.log"
    params:
        header="Molecule\tCount"
    shell:
        "samtools idxstats {input.bam} | awk 'BEGIN{{OFS='\\t'}} $3 != 0' - | cut -f1,3 -| sed '1 i\{params.header}' - > {output.counts}"

# Count deduplicated read counts from BAM file and extract relevant columns from output.
# Results will be written to tab-delimited .tsv file.
rule count_deduplicated_read_counts:
    input:
        bam="results/{sample}/bam/{sample}_merged_dedup.bam",
        bai="results/{sample}/bam/{sample}_merged_dedup.bai"
    output:
        counts="results/{sample}/{sample}_counts_dedup.tsv"
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/count_deduplicated_read_counts.log"
    shell:
        "samtools idxstats {input.bam} | awk 'BEGIN{{OFS='\\t'}} $3 != 0' - | cut -f1,3 -| sed '1 i\Molecule\tCount' - > {output.counts}"

# Merge sample-specific raw read counts to one, merged .tsv file where each 
# column represents one sample. External Python script will be used for merging.
rule merge_raw_read_counts:
    input:
        expand("results/{sample}/{sample}_counts.tsv", sample=get_wildcards())
    output:
        "results/counts.tsv"
    conda:
        "../envs/env.yaml"
    log:
        expand("results/{sample}/logs/merge_raw_read_counts.log", sample=get_wildcards())
    script:
        "../scripts/merge_read_counts.py"

# Merge sample-specific and deduplicated read counts to one, merged .tsv file 
# where each column represents one sample. External Python script will be used for merging.
rule merge_deduplicated_read_counts:
    input:
        expand("results/{sample}/{sample}_counts_dedup.tsv", sample=get_wildcards())
    output:
        "results/counts_deduplicated.tsv"
    conda:
        "../envs/env.yaml"
    log:
        expand("results/{sample}/logs/merge_deduplicated_read_counts.log", sample=get_wildcards())
    script:
        "../scripts/merge_read_counts.py"