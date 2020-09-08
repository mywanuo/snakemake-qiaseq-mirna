# Map reads to mirBase sequences and convert directly to sorted .bam file.
# Non-aligned reads will be saved as a new .fastq file and be forwarded to 
# alignment step (piRNA alignment). Index file (.bai) for sorted reads
# will be generated after the alignment.
rule map_mirna_sequences:
    input:
        fastq="results/{sample}/fastq/{sample}_trimmed.fastq.gz",
        fasta="resources/mirna.fasta",
        index=[
            "resources/mirna.1.bt2",
            "resources/mirna.2.bt2",
            "resources/mirna.3.bt2",
            "resources/mirna.4.bt2",
            "resources/mirna.rev.1.bt2",
            "resources/mirna.rev.2.bt2"
        ]
    output:
        bam="results/{sample}/bam/{sample}_mirna.bam",
        bai="results/{sample}/bam/{sample}_mirna.bai",
        fastq=temp("results/{sample}/fastq/{sample}_mirna_unmapped.fastq.gz")
    threads:
        4
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/map_mirna_sequences.log"
    params:
        index="resources/mirna"
    threads:
        4
    shell:
        """
        bowtie2 --very-sensitive-local -p {threads} -x {params.index} --un-gz {output.fastq} {input.fastq} | samtools sort - > {output.bam} 2> {log} 
        samtools index {output.bam} {output.bai}
        """

# Map reads to piRNAdb sequences and convert directly to sorted .bam file.
# Non-aligned reads will be saved as a new .fastq file and be forwarded to 
# alignment step (reference genome alignment). Index file (.bai) for sorted 
# reads will be generated after the alignment.
rule map_pirna_sequences:
    input:
        fastq="results/{sample}/fastq/{sample}_mirna_unmapped.fastq.gz",
        fasta="resources/pirna.fasta",
        index=[
            "resources/pirna.1.bt2",
            "resources/pirna.2.bt2",
            "resources/pirna.3.bt2",
            "resources/pirna.4.bt2",
            "resources/pirna.rev.1.bt2",
            "resources/pirna.rev.2.bt2"
        ]
    output:
        bam="results/{sample}/bam/{sample}_pirna.bam",
        bai="results/{sample}/bam/{sample}_pirna.bai",
        fastq="results/{sample}/fastq/{sample}_unmapped.fastq.gz"
    threads:
        4
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/map_pirna_sequences.log"
    params:
        index="resources/pirna"
    threads:
        4
    shell:
        """
        bowtie2 --very-sensitive-local -p {threads} -x {params.index} --un-gz {output.fastq} {input.fastq} | samtools sort - > {output.bam} 2> {log} 
        samtools index {output.bam} {output.bai}
        """

# Merge individual BAM files to one, merged BAM file to reduce the number of 
# processed files and simplify following workflow.
rule merge_bam_files:
    input:
        bam=[
            "results/{sample}/bam/{sample}_mirna.bam",
            "results/{sample}/bam/{sample}_pirna.bam"
        ],
        bai=[
            "results/{sample}/bam/{sample}_mirna.bai",
            "results/{sample}/bam/{sample}_pirna.bai"
        ] 
    output:
        bam=temp("results/{sample}/bam/{sample}_merged.bam"),
        bai=temp("results/{sample}/bam/{sample}_merged.bai")
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/merge_bam_files.log"
    threads:
        4
    shell:
        """
        samtools merge --threads {threads} {output.bam} {input.bam} > {log} 
        samtools index {output.bam} {output.bai}  
        """