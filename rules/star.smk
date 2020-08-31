rule star_index:
    input:
        fasta = get_fasta_path(),
        gtf = get_gtf_path(),
        index = get_index_path()
    output:
        directory("star/index/idx")
    message:
        "Indexing idx with STAR"
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(35840 + 10240 * attempt, 51200)
        ),
        time_min = (
            lambda wildcards, attempt: min(35 * attempt, 120)
        )
    params:
        gtf = "resources/idx.gtf",
        sjdbOverhang = "100",
        extra = ""
    log:
        "logs/star/index.log"
    wrapper:
        f"{git}/bio/star/index"


rule star_mapping:
    input:
        unpack(star_sample_pair_w),
        index = "star/index/idx"
    output:
        temp("star/bam/{sample}/Aligned.sortedByCoord.out.bam")
    message:
        "Mapping {wildcards.sample} with STAR"
    wildcard_constraints:
        sample = "|".join(design.Sample_id)
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(35840 + 10240 * attempt, 51200)
        ),
        time_min = (
            lambda wildcards, attempt: min(50 * attempt, 120)
        )
    params:
        index = "star/index/idx",
        extra = ("--outSAMtype BAM SortedByCoordinate "
                 "--outSAMattributes All "
                 "--outFilterType BySJout "
                 "--outFilterMismatchNmax 999 "
                 "--alignSJDBoverhangMin 1 "
                 "--outFilterMismatchNoverReadLmax 0.04 "
                 "--alignMatesGapMax 1000000 "
                 "--alignIntronMax 1000000 "
                 "--alignIntronMin 20 "
                 "--alignSJoverhangMin 8 "
                 "--outFilterMultimapNmax 20 "
                 " --quantMode GeneCounts ")
    log:
        "logs/star/bam/{sample}.log"
    wrapper:
        f"{git}/bio/star/align"


rule star_rename:
    input:
        "star/bam/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "star/bam/{sample}.bam"
    message:
        "Renaming {wildcards.sample} for further analyses"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(5 * attempt, 15)
        )
    params:
        "--verbose"
    conda:
        "../envs/bash.yaml"
    log:
        "logs/rename/{sample}.log"
    shell:
        "cp {params} {input} {output} > {log} 2>&1"


rule samtools_bam_index:
    input:
        "star/bam/{sample}.bam"
    output:
        "star/bam/{sample}.bam.bai"
    message:
        "Indexing {wildcards.sample} bam file"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(512 * attempt, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(15 * attempt, 45)
        )
    params:
        ""
    wrapper:
        f"{git}/bio/samtools/index"


rule samtools_flagstat:
    input:
        "star/bam/{sample}.bam"
    output:
        temp("samtools/flagstat/{sample}.flagstat")
    message:
        "Gathering statistics on {wildcards.sample}'s mapping"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(512 * attempt, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(15 * attempt, 45)
        )
    params:
        ""
    wrapper:
        f"{git}/bio/samtools/flagstat"
