rule all:
	input:
		expand('{sample}_R1.fastq.gz', sample=SAMPLES),
		expand('{sample}_R1_temp1.fastq.gz', sample=SAMPLES),
		expand('{sample}_R1_temp2.fastq.gz', sample=SAMPLES),
		expand('{sample}_R1_temp3.fastq.gz', sample=SAMPLES),
		expand('{sample}_R1_trimmed.fastq.gz', sample=SAMPLES),
		expand('{sample}_mirna.bam', sample=SAMPLES),
		expand('{sample}_mirna.bai', sample=SAMPLES),
		expand('{sample}_mirna_dedup.bam', sample=SAMPLES),
		expand('{sample}_mirna_dedup.bai', sample=SAMPLES),
		expand('{sample}_mirna_counts.tsv', sample=SAMPLES),
        expand('merged.tsv')