#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script aims to prepare the list of files to be processed
by the infer-sequencing-libraries pipeline

It iterates over a given directory, lists all fastq files. As a pair of
fastq files usually have names that follows each other in the alphabetical
order, this script sorts the fastq files names and, by default, creates
pairs of fastq files that way.

Finally, it writes these pairs, using the longest common substring as
identifier. The written file is a TSV file.

Reference files' paths are stored in a yaml-formatted configuration file.

You can test this script with:
pytest -v ./prepare_pipeline.py

Usage example:
# Single ended reads example:
python3.8 ./prepare_pipeline.py tests/reads tests/genome/genome.chr21.fa tests/genome/annotation.chr21.gtf tests/genome/chr21.bed --single

# Paired-end libary example:
python3.8 ./prepare_pipeline.py tests/reads tests/genome/genome.chr21.fa tests/genome/annotation.chr21.gtf tests/genome/chr21.bed
"""

import argparse  # Parse command line
import logging  # Traces and loggings
import os  # OS related activities
import pandas as pd  # Parse TSV files
import pytest  # Unit testing
import shlex  # Lexical analysis
import sys  # System related methods
import yaml  # Handle YAML formatting

from pathlib import Path  # Paths related methods
from snakemake.utils import makedirs  # Easily build directories
from typing import Dict, Generator, List, Any  # Type hints

from common_script_infer_sequencing_libraries import CustomFormatter


# Processing functions
# Looking for fastq files
def search_fq(
    fq_dir: Path, recursive: bool = False
) -> Generator[str, str, None]:
    """
    Iterate over a directory and search for fastq files

    Parameters:
        fq_dir      Path        Path to the fastq directory in which to search
        recursive   bool        A boolean, weather to search recursively in
                                sub-directories (True) or not (False)

    Return:
                    Generator[str, str, None]       A Generator of paths

    Example:
    >>> search_fq(Path("tests/reads/"))
    <generator object search_fq at 0xXXXXXXXXXXXX>

    >>> list(search_fq(Path("tests/", True)))
    [PosixPath('tests/reads/a.chr21.2.fastq'),
     PosixPath('tests/reads/b.chr21.2.fastq'),
     PosixPath('tests/reads/a.chr21.1.fastq'),
     PosixPath('tests/reads/b.chr21.1.fastq')]
    """
    for path in fq_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_fq(path, recursive)
            else:
                continue

        if path.name.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
            yield path


# Testing search_fq
def test_search_fq():
    """
    This function tests the ability of the function "search_fq" to find the
    fastq files in the given directory

    Example:
    pytest -v prepare_design.py -k test_search_fq
    """
    path = Path("tests/reads/")

    expected = list(
        path / "{}.chr21.{}.fq.gz".format(sample, stream)
        for sample in ["a", "b"]
        for stream in [1, 2]
    )
    assert sorted(list(search_fq(path))) == sorted(expected)


# Turning the FQ list into a dictionnary
def classify_fq(fq_files: List[Path], paired: bool = True) -> Dict[str, Path]:
    """
    Return a dictionnary with identified fastq files (paried or not)

    Parameters:
        fq_files    List[Path]      A list of paths to iterate over
        paired      bool            A boolean, weather the dataset is
                                    pair-ended (True) or single-ended (False)

    Return:
                    Dict[str, Path] A dictionnary: for each Sample ID, the ID
                                    is repeated alongside with the upstream
                                    /downstream fastq files.

    Example:
    # Paired-end single sample
    >>> classify_fq([Path("file1.R1.fq"), Path("file1.R2.fq")], True)
    {'file1.R1.fq': {'Downstream_file': PosixPath("/path/to/file1.R1.fq"),
     'Sample_id': 'file1.R1',
     'Upstream_file': PosixPath('/path/to/file1.R2.fq')}

    # Single-ended single sample
    >>> classify_fq([Path("file1.fq")], False)
    {'file1.fq': {'Sample_id': 'file1',
     'Upstream_file': PosixPath('/path/to/file1.fq')}}
    """
    fq_dict = {}
    if paired is not True:
        # Case single fastq per sample
        logging.debug("Sorting fastq files as single-ended")
        for fq in fq_files:
            fq_dict[fq.name] = {
                "Sample_id": fq.stem,
                "Upstream_file": fq.absolute(),
            }
    else:
        # Case pairs of fastq are used
        logging.debug("Sorting fastq files as pair-ended")
        for fq1, fq2 in zip(fq_files[0::2], fq_files[1::2]):
            fq_dict[fq1.name] = {
                "Sample_id": fq1.stem,
                "Upstream_file": fq1.absolute(),
                "Downstream_file": fq2.absolute(),
            }
    logging.debug(fq_dict)
    return fq_dict


def test_classify_fq():
    """
    This function takes input from the pytest decorator
    to test the classify_fq function

    Example:
    pytest -v ./prepare_design.py -k test_classify_fq
    """
    prefix = Path(__file__).parent.parent
    fq_list = sorted(list(search_fq(prefix / "tests" / "reads")))
    expected = {
        "a.chr21.1.fq.gz": {
            "Sample_id": "a.chr21.1.fq",
            "Upstream_file": prefix / "tests" / "reads" / "a.chr21.1.fq.gz",
            "Downstream_file": prefix / "tests" / "reads" / "a.chr21.2.fq.gz",
        },
        "b.chr21.1.fq.gz": {
            "Sample_id": "b.chr21.1.fq",
            "Upstream_file": prefix / "tests" / "reads" / "b.chr21.1.fq.gz",
            "Downstream_file": prefix / "tests" / "reads" / "b.chr21.2.fq.gz",
        },
    }

    assert classify_fq(fq_list) == expected


# Parsing command line arguments
# This function won't be tested
def parse_args(args: Any = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    Build a command line parser object

    Parameters:
        args    Any                 Command line arguments

    Return:
                ArgumentParser      Parsed command line object

    Example:
    >>> parse_args(shlex.split("/path/to/fasta --single"))
    Namespace(debug=False, output='design.tsv', path='/path/to/fasta',
    quiet=False, recursive=False, single=True)
    """
    # Defining command line options
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does not perform any magic. Check the result.",
    )

    # Required arguments
    main_parser.add_argument(
        "path",
        help="Path to the directory containing fastq files",
        type=str
    )

    main_parser.add_argument(
        "fasta",
        help="Path to fasta-formatted genome sequence",
        type=str
    )

    main_parser.add_argument(
        "gtf",
        help="Path to GTF-formatted genome annotation",
        type=str
    )

    main_parser.add_argument(
        "bed",
        help="Path to BED12-formatted genome reference (see RSeQC)",
        type=str
    )

    # Optional arguments
    main_parser.add_argument(
        "-s",
        "--single",
        help="The samples are single ended rnaseq reads, not pair ended",
        action="store_true",
    )

    main_parser.add_argument(
        "-r",
        "--recursive",
        help="Recursively search in sub-directories for fastq files",
        action="store_true",
    )

    main_parser.add_argument(
        "--design",
        help="Path to output design file (default: %(default)s)",
        type=str,
        default="design.tsv",
    )

    main_parser.add_argument(
        "--config",
        help="Path to output config file (default: %(default)s)",
        type=str,
        default="config.yaml",
    )

    # Logging options
    log = main_parser.add_mutually_exclusive_group()
    log.add_argument(
        "-d",
        "--debug",
        help="Set logging in debug mode",
        default=False,
        action="store_true",
    )

    log.add_argument(
        "-q",
        "--quiet",
        help="Turn off logging behaviour",
        default=False,
        action="store_true",
    )

    # Parsing command lines
    return main_parser.parse_args(args)


def test_parse_args() -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    options = parse_args(shlex.split(
        "/path/to/fastq/dir/ "
        "/path/to/fasta "
        "/path/to/gtf "
        "/path/to/bed"
    ))
    expected = argparse.Namespace(
        debug=False,
        design="design.tsv",
        config="config.yaml",
        path="/path/to/fastq/dir/",
        fasta="/path/to/fasta",
        gtf="/path/to/gtf",
        bed="/path/to/bed",
        quiet=False,
        recursive=False,
        single=False,
    )
    assert options == expected


# Main function, the core of this script
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole preparation sequence

    Parameters:
        args    ArgumentParser      The parsed command line

    Example:
    >>> main(parse_args(shlex.split("/path/to/fasta/dir/")))
    """
    # Searching for fastq files and sorting them alphabetically
    fq_files = sorted(list(search_fq(Path(args.path), args.recursive)))
    logging.debug("Head of alphabeticaly sorted list of fastq files:")
    logging.debug([str(i) for i in fq_files[0:5]])

    # Building a dictionnary of fastq (pairs?) and identifiers
    fq_dict = classify_fq(fq_files)

    # Using Pandas to handle TSV output (yes pretty harsh I know)
    data = pd.DataFrame(fq_dict).T
    logging.debug("\n{}".format(data.head()))
    logging.debug("Saving design to {}".format(args.design))
    data.to_csv(args.design, sep="\t", index=False)

    config = {
        "fasta": args.fasta,
        "gtf": args.gtf,
        "bed": args.bed
    }
    logging.debug("Saving configuration to {}".format(args.config))
    logging.debug("Configurations: {}".format(str(config)))
    with open(args.config, "w") as yamlout:
        yamlout.write(yaml.dump(config))


# Running programm if not imported
if __name__ == "__main__":
    # Parsing command line
    args = parse_args()

    makedirs("logs/prepare")

    # Build logging object and behaviour
    logging.basicConfig(
        filename="logs/prepare/design.log", filemode="w", level=logging.DEBUG
    )

    try:
        logging.debug("Preparing design")
        main(args)
    except Exception as e:
        logging.exception("%s", e)
        raise
    sys.exit(0)
