import sys
from snakemake.utils import min_version

if sys.version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later")

min_version('5.22.1')
container: "docker://continuumio/miniconda3:5.0.1"

include: "rules/common.smk"
include: "rules/copy.smk"
include: "rules/fastq.smk"
include: "rules/picard.smk"
include: "rules/seqt.smk"
include: "rules/star.smk"
include: "rules/rseqc.smk"

rule all:
    input:
        **get_targets(get_manuals=True, get_stats=True, get_rcqc=True)
    message:
        "Finishing pipeline"
