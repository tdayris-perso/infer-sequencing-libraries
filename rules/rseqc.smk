rule infer_experiment:
    input:
        bed = get_bed_genome_path()
        bam = "star/bam/{sample}.bam"
        bam_index = "star/bam/{sample}.bam.bai"
    output:
        "rseqc/infer_experiment/{sample}.txt"
    message:
        "Gathering experiment inference from RSeQC on {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2 * 128,
        time_min = lambda wildcards, attempt: attempt * 15
    log:
        "logs/rseqc/infer_experiment/{sample}.log"
    shell:
        " infer_experiment.py "
        " -r {input.bed} "
        " -i {input.bam} "
        " > {output} "
        " 2> {log} "
