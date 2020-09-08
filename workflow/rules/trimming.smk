# Merge separate fastq files into one file.
rule merge_fastq_files:
    input:
        get_fastq
    output:
        "results/{sample}/fastq/{sample}_untrimmed.fastq.gz"
    shell:
        "cat {input} > {output}" 

# Trim Illumina Universal Adapter sequences from the 3' end of reads.   
rule trim_illumina_adapters:
    input:
        "results/{sample}/fastq/{sample}_untrimmed.fastq.gz"
    output:
        temp("results/{sample}/fastq/{sample}_illumina_trimmed.fastq.gz")
    params:
        adapter=config["trimming"]["illumina_adapter"]
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/trim_illumina_adapters.log"
    threads:
        4
    shell:
        "cutadapt -a {params.adapter} --discard-untrimmed -m 20 -q 10 -j {threads} -o {output} --quiet {input} > {log}"

# Trim UMI sequences from the 3' end of trimmed read and add them to read header.
rule trim_umi_sequences:
    input:
        "results/{sample}/fastq/{sample}_illumina_trimmed.fastq.gz"
    output:
        temp("results/{sample}/fastq/{sample}_umi_trimmed.fastq.gz")
    params:
        umi_sequence=config["trimming"]["umi_sequence"]
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/trim_umi_sequences.log"
    shell:
        "umi_tools extract --3prime --stdin={input} --bc-pattern={params.umi_sequence} --stdout={output} > {log} "

# Trim 5' adapter sequences from the 5' end of the read.
rule trim_5_adapter:
    input:
        "results/{sample}/fastq/{sample}_umi_trimmed.fastq.gz"
    output:
        temp("results/{sample}/fastq/{sample}_5adapter_trimmed.fastq.gz")
    params:
        adapter=config["trimming"]["adapter_5"]
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/trim_5_adapter.log"
    threads:
        4
    shell:
        "cutadapt -g {params.adapter} -m 20 -j {threads} -o {output} --quiet {input} > {log}"

# Trim 3' adapter sequence from the 3' end of the read.
# Discard untrimmed reads and perform some quality-trimming.
rule trim_3_adapter:
    input:
        "results/{sample}/fastq/{sample}_5adapter_trimmed.fastq.gz"
    output:
        "results/{sample}/fastq/{sample}_trimmed.fastq.gz"
    params:
        adapter_3=config["trimming"]["adapter_3"],
        min_length=config["trimming"]["min_read_length"],
        min_quality=config["trimming"]["min_read_quality"] 
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/trim_3_adapter.log"
    threads:
        4
    shell:
        "cutadapt -a {params.adapter_3} --discard-untrimmed -m {params.min_length} -q {params.min_quality} -j {threads} -o {output} --quiet {input} > {log}"
