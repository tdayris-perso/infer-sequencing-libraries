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
