rule read_length:
    input:
        expand(
            "raw_data/{sample}.fastq.gz",
            stream = ["1", "2"],
            allow_missing = True
        )
    output:
        temp("read_length/{sample}.txt")
    message:
        "Gathering read length information on {wildcards.sample}"
    group:
        "Manual_data_gathering"
    threads:
        1
    resources:
        time_min = lambda wildcards, attempt: attempt * 30,
        mem_mb = lambda wildcards, attempt: 256
    log:
        "logs/read/length/{sample}.log"
    script:
        "scripts/read_length.py"


rule read_quality:
    input:
        "raw_data/{sample}.fastq.gz"
    output:
        temp("read_quality/{sample}.txt")
    message:
        "Gathering read qualities information on {wildcards.sample}"
    group:
        "Manual_data_gathering"
    threads:
        1
    resources:
        time_min = lambda wildcards, attempt: attempt * 30,
        mem_mb = lambda wildcards, attempt: 256
    log:
        "logs/read/quality/{sample}.log"
    script:
        "scripts/read_quality.py"


rule read_number:
    input:
        "raw_data/{sample}.fastq.gz"
    output:
        temp("read_number/{sample}.txt")
    message:
        "Gathering read number information on {wildcards.sample}"
    group:
        "Manual_data_gathering"
    threads:
        1
    resources:
        time_min = 15,
        mem_mb = 128
    log:
        "logs/read/number/{sample}.log"
    shell:
        " awk 'END {{ print NR/4 }}' {input} > {output} 2> {log}"


rule paste_read_information:
    input:
        expand(
            "read_{data}/{sample}.txt",
            data=["length", "number", "quality"]
            allow_missing=True
        )
    output:
        temp("stats/manual/{sample}.txt")
    message:
        "Gathering all metrics for {wildcards.sample}"
    group:
        "Manual_data_gathering"
    threads:
        1
    resources:
        time_min = 5,
        mem_mb = 128
    log:
        "logs/read/sats/{sample}.log"
    shell:
        "paste {input} > {output} 2> {log}"


rule cat_read_information:
    input:
        expand(
            "stats/manual/{sample}.txt",
            sample=design.Sample_id
        )
    output:
        "stats/manual/complete.txt"
    message:
        "Gathering information for all samples"
    threads:
        1
    resources:
        time_min = lambda wildcards, attempt: attempt * 30,
        mem_mb = lambda wildcards, attempt: 256
    log:
        "logs/read/complete.log"
    shell:
        "cat {input} > {output} 2> {log}"
