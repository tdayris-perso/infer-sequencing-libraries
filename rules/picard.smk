rule create_dict:
    input:
        get_fasta_path()
    output:
        get_dict_path()
    log:
        "logs/picard/create_dict.log"
    params:
        extra=""  # optional: extra arguments for picard.
    wrapper:
        f"{git}/bio/picard/createsequencedictionary"


rule samtools_index:
    input:
        get_fasta_path()
    output:
        get_index_path()
    params:
        "" # optional params string
    wrapper:
        f"{git}/bio/samtools/faidx"


rule alignment_summary:
    input:
        ref=get_fasta_path(),
        ref_idx=get_index_path(),
        ref_dict=get_dict_path(),
        bam="star/bam/{sample}.bam"
    output:
        "picard/collectalignmentsummarymetrics/{sample}.summary.txt"
    log:
        "logs/picard/alignment-summary/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        f"{git}/bio/picard/collectalignmentsummarymetrics"


rule insert_size:
    input:
        "star/bam/{sample}.bam"
    output:
        txt="picard/collectinsertsizemetrics/{sample}.isize.txt",
        pdf="picard/collectinsertsizemetrics/{sample}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        f"{git}/bio/picard/collectinsertsizemetrics"


rule cat_alignment_summary_metrix:
    input:
        "picard/collectalignmentsummarymetrics/{sample}.summary.txt"
    output:
        temp("picard/collectalignmentsummarymetrics/{sample}.processed.tsv")
    message:
        "Processing alignment summary for {wildcards.sample}"
    threads:
        1
    resources:
        time_min = 10,
        mem_mb = 512
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/pandas/collectalignmentsummarymetrics/format_{sample}.log"
    script:
        "../scripts/picard_summary_to_tsv.py"


rule read_alignment_summary_metrix:
    input:
        expand(
            "picard/collectalignmentsummarymetrics/{sample}.processed.tsv",
            sample=design.Sample_id
        )
    output:
        temp("stats/collectalignmentsummarymetrics/complete.tsv")
    message:
        "Collecting all summaries"
    threads:
        1
    resources:
        time_min = 10,
        mem_mb = 512
    conda:
        "../envs/bash.yaml"
    log:
        "logs/pandas/collectalignmentsummarymetrics/cat.log"
    shell:
        ' echo -e '
        '"Sample_id\tUpstream_reads\tDownstream_reads\tLibrary" '
        ' > {output} 2> {log} && '
        ' cat {input} >> {output} 2>> {log} '


rule cat_picard_insert_size_metrics:
    input:
        "picard/collectinsertsizemetrics/{sample}.isize.txt"
    output:
        temp("picard/collectinsertsizemetrics/{sample}.isize.tsv")
    message:
        "Processing insert sizes for {wildcards.sample}"
    threads:
        1
    resources:
        time_min = 10,
        mem_mb = 512
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/pandas/collectinsertsizemetrics/format_{sample}.log"
    script:
        "../scripts/picard_insert_size_to_tsv.py"


rule read_picard_insert_size_metrics:
    input:
        expand(
            "picard/collectinsertsizemetrics/{sample}.isize.tsv",
            sample=design.Sample_id
        )
    output:
        temp("stats/collectinsertsizemetrics/complete.isize.tsv")
    message:
        "Gathering all insert sizes"
    threads:
        1
    resources:
        time_min = 10,
        mem_mb = 512
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/pandas/collectinsertsizemetrics/cat.log"
    shell:
        ' echo -e '
        ' "Sample_id\tMean_insert_size\tStd_insert_size\tOrientation" '
        ' > {output} 2> {log} && '
        ' cat {input} >> {output} 2>> {log} '
