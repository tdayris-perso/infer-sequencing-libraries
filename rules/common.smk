"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""

import os.path
import pandas

from typing import Any, List           # Type hinting
from snakemake.utils import validate   # Check Yaml/TSV formats

from common_ngs_cleaning import sample_stream, fastq_pairs

# github prefix
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/Unofficial"


# Loading configuration
if config == dict():
    configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pandas.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")


ref_link_dict = {
    os.path.basename(config["ref"]["fasta"]): config["ref"]["fasta"],
    os.path.basename(config["ref"]["gtf"]): config["ref"]["gtf"],
    os.path.basename(config["ref"]["bed"]): config["ref"]["bed"]
}


def star_sample_pair_w(wildcards: Any) -> Dict[str, str]:
    return {
        "fq1": f"seqt/{wildcards.sample}.1.fastq.gz",
        "fq2": f"seqt/{wildcards.sample}.2.fastq.gz"
    }


def get_fasta_path() -> str:
    """
    Return copied fasta path
    """
    base = os.path.basename(config["ref"]["fasta"])
    return f"genomes/{base}"


def get_index_path() -> str:
    """
    Return copied fasta index path
    """
    base = os.path.basename(config["ref"]["fasta"])
    return f"genomes/{base}.fai"


def get_dict_path() -> str:
    """
    Return copied fasta sequence dictionary path
    """
    base = os.path.basename(config["ref"]["fasta"])
    noext = ".".join(os.path.splitext(base)[:-1])
    return f"genomes/{noext}.dict"


def get_bed_genome_path() -> str:
    """
    Return copied bed12 genome
    """
    base = os.path.basename(config["ref"]["bed"])
    return f"genomes/{base}"


def get_gtf_path() -> str:
    """
    Return copied gtf path
    """
    base = os.path.basename(config["ref"]["gtf"])
    return f"genomes/{base}"
