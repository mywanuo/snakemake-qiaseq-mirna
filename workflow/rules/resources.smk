rule get_mirna_sequences:
    output:
        mature=temp("resources/mirna_mature.fasta"),
        hairpin=temp("resources/mirna_hairpin.fasta"),
        merged="resources/mirna.fasta"
    conda:
        "../envs/env.yaml"
    params:
        url_mature=config["fasta"]["mature_mirna"],
        url_hairpin=config["fasta"]["hairpin_mirna"] 
    log:
        "resources/logs/get_mirna_sequences.log"
    shell:
        """
        curl -s {params.url_mature} | zcat - | seqkit grep -r -p ^hsa | seqkit seq --rna2dna > {output.mature} 2> {log} 
        curl -s {params.url_hairpin} | zcat - | seqkit grep -r -p ^hsa | seqkit seq --rna2dna > {output.hairpin} 2> {log} 
        cat {output.mature} {output.hairpin}  > {output.merged}
        """

rule get_mirna_index:
    input:
        fasta="resources/mirna.fasta"
    output:
        index=[
            "resources/mirna.1.bt2",
            "resources/mirna.2.bt2",
            "resources/mirna.3.bt2",
            "resources/mirna.4.bt2",
            "resources/mirna.rev.1.bt2",
            "resources/mirna.rev.2.bt2"
        ]
    conda:
        "../envs/env.yaml"
    log:
        "resources/logs/get_mirna_index.log"
    params:
        prefix="resources/mirna"
    threads:
        4
    shell:
        """
        bowtie2-build --threads {threads} --quiet {input.fasta} {params.prefix} > {log} 
        """

rule get_pirna_sequences:
    output:
        "resources/pirna.fasta"
    conda:
        "../envs/env.yaml"
    log:
        "resources/logs/get_pirna_index.log"
    params:
        url=config["fasta"]["pirna"]
    shell:
        "curl -s {params.url} | zcat - > {output} 2> {log}"

rule get_pirna_index:
    input:
        "resources/pirna.fasta"
    output:
        "resources/pirna.1.bt2",
        "resources/pirna.2.bt2",
        "resources/pirna.3.bt2",
        "resources/pirna.4.bt2",
        "resources/pirna.rev.1.bt2",
        "resources/pirna.rev.2.bt2"
    conda:
        "../envs/env.yaml"
    log:
        "resources/logs/get_pirna_index.log"
    params:
        prefix="resources/pirna"
    threads:
        4
    shell:
        "bowtie2-build --threads {threads} --quiet {input} {params.prefix} > {log}"