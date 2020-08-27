rule copy_extra:
    input:
        lambda wildcards: ref_link_dict[wildcards.files]
    output:
        temp("genomes/{files}")
    message:
        "Copying {wildcards.files} as reference"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 128, 512),
        time_min = lambda wildcards, attempt: min(attempt * 1440, 2832)
    log:
        "logs/copy/{files}.log"
    wildcard_constraints:
        files = r"[^/]+"
    threads: 1
    params:
        extra = config["params"].get("copy_extra", ""),
        cold_storage = config.get("cold_storage", ["NONE"])
    wrapper:
        f"{git}/bio/cp"