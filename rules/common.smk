"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""

import os.path
import pandas

from typing import Any, Dict, Union, List  # Type hinting
from snakemake.utils import validate  # Check Yaml/TSV formats

# github prefix
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/Unofficial"


# Loading configuration
if config == dict():
    configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pandas.read_csv(
    "design.tsv",
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")


ref_link_dict = {
    os.path.basename(config["fasta"]): config["fasta"],
    os.path.basename(config["gtf"]): config["gtf"],
    os.path.basename(config["bed"]): config["bed"]
}

is_paired = False
fq_link_dict = {}
if "Downstream_file" in design.columns.tolist():
    is_paired = True
    zip_fq = zip(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )
    for sample_id, up, down in zip_fq:
        fq_link_dict[f"{sample_id}_R1"] = up
        fq_link_dict[f"{sample_id}_R2"] = down
else:
    zip_fq = zip(
        design.Sample_id,
        design.Upstream_file
    )
    for sample_id, up, down in zip_fq:
        fq_link_dict[f"{sample_id}"] = up


def get_samples_w(wildcards: Any) -> List[str]:
    """
    Return a single file for single-ended libraries
    or pairs for paired-end ones
    """
    if is_paired is True:
        return [
            fq_link_dict[f"{wildcards.sample}_R1"],
            fq_link_dict[f"{wildcards.sample}_R2"]
        ]
    return [
        fq_link_dict[wildcards.sample]
    ]


def get_seqt_pairs_w(wildcards: Any) -> Union[str, Dict[str, str]]:
    """
    Return fastq files for seqtk
    """
    if is_paired is True:
        return {
            "f1": fq_link_dict[f"{wildcards.sample}_R1"],
            "f2": fq_link_dict[f"{wildcards.sample}_R2"]
        }
    return fq_link_dict[wildcards.sample]


def star_sample_pair_w(wildcards: Any) -> Dict[str, str]:
    return {
        "fq1": f"seqt/{wildcards.sample}.1.fastq.gz",
        "fq2": f"seqt/{wildcards.sample}.2.fastq.gz"
    }


def get_fasta_path() -> str:
    """
    Return copied fasta path
    """
    base = os.path.basename(config["fasta"])
    return f"genomes/{base}"


def get_index_path() -> str:
    """
    Return copied fasta index path
    """
    base = os.path.basename(config["fasta"])
    return f"genomes/{base}.fai"


def get_dict_path() -> str:
    """
    Return copied fasta sequence dictionary path
    """
    base = os.path.basename(config["fasta"])
    noext = ".".join(os.path.splitext(base)[:-1])
    return f"genomes/{noext}.dict"


def get_bed_genome_path() -> str:
    """
    Return copied bed12 genome
    """
    base = os.path.basename(config["bed"])
    return f"genomes/{base}"


def get_gtf_path() -> str:
    """
    Return copied gtf path
    """
    base = os.path.basename(config["gtf"])
    return f"genomes/{base}"


def get_targets(get_manuals: bool = False,
                get_bam: bool = False) -> Dict[str, str]:
    """
    Return all needed output files
    """
    targets = {}

    if get_manuals is True:
        targets["manuals"] = "stats/manual/complete.txt"


    if get_bam is True:
        targets["bam"] = expand(
            "star/bam/{sample}.bam",
            sample = design.Sample_id
        )

    return targets
