rule subset_reads:
    input:
        unpack(get_seqt_pairs_w)
    output:
        f1 = "seqt/{sample}.1.fastq.gz",
        f2 = "seqt/{sample}.2.fastq.gz"
    message:
        "Subsampling {wildcards.sample}"
    threads:
        2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4 * 1024,
        time_min = lambda wildcards, attempt: attempt * 30
    params:
        n = config["params"].get("seqt_read_number", 100000),
        seed = config["params"].get("seqt_seed", 1234)
    log:
        "logs/seqt/{sample}.log"
    wrapper:
        f"{git}/seqtk/subsample/pe"
