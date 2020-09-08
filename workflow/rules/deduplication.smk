# Deduplicate reads on the basis of previously identified UMI sequences with 
# UMItools dedup command.
rule deduplicate_reads:
    input:
        bam="results/{sample}/bam/{sample}_merged.bam",
        bai="results/{sample}/bam/{sample}_merged.bai"
    output:
        bam="results/{sample}/bam/{sample}_merged_dedup.bam",
        bai="results/{sample}/bam/{sample}_merged_dedup.bai"
    conda:
        "../envs/env.yaml"
    log:
        "results/{sample}/logs/deduplicate_reads.log"
    shell:
        """
        umi_tools dedup -I {input.bam} -S {output.bam} > {log} 
        samtools index {output.bam} {output.bai}
        """