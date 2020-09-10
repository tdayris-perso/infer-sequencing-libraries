rule rcqc_infer_experiment:
    input:
        bed = get_bed_genome_path(),
        bam = "star/bam/{sample}.bam",
        bam_index = "star/bam/{sample}.bam.bai"
    output:
        "rseqc/infer_experiment/{sample}.txt"
    message:
        "Gathering experiment inference from RSeQC on {wildcards.sample}"
    group:
        "RSeQC"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2 * 128,
        time_min = lambda wildcards, attempt: attempt * 15
    log:
        "logs/rseqc/infer_experiment/{sample}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        " infer_experiment.py "
        " -r {input.bed} "
        " -i {input.bam} "
        " > {output} "
        " 2> {log} "


rule rseqc_format_result:
    input:
        "rseqc/infer_experiment/{sample}.txt"
    output:
        "rseqc/infer_experiment/{sample}.tsv"
    message:
        "Formatting rseqc output"
    group:
        "RSeQC"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2 * 128,
        time_min = lambda wildcards, attempt: attempt * 15
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/rseqc/format/{sample}.log"
    script:
        "../scripts/rseqc_to_tsv.py"


rule rseqc_merge_results:
    input:
        expand(
            "rseqc/infer_experiment/{sample}.tsv",
            sample=design.Sample_id
        )
    output:
        "rseqc/merged.tsv"
    message:
        "Concatenating all rseqc results"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2 * 128,
        time_min = lambda wildcards, attempt: attempt * 15
    conda:
        "../envs/bash.yaml"
    log:
        "log/rseqc/merge.log"
    shell:
        "for i in {input}; do sed '1d' {input}; done > {output} 2> {log}"
