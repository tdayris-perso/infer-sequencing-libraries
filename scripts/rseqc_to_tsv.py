#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script parses RSeQC results
"""

import logging

logging.basicConfig(
    filename=snakemake.log[0], filemode="w", level=logging.DEBUG
)

titles = "\t".join([
    "Sample_id",
    "Undetermined_Fraction",
    "FR_Fraction",
    "RF_Fraction",
])

with open(snakemake.input[0]) as rseqc:
    fractions = [
        line.split(":")[-1][:-1]
        for line in rseqc
        if not (line == "\n" or line.startswith("This is"))
    ]
logging.debug(f"{snakemake.wildcards.sample} processed")

content = "\t".join([snakemake.wildcards.sample] + fractions)

with open(snakemake.output[0], 'w') as outtsv:
    outtsv.write("\n".join([titles, content]))
    logging.debug("Process over")
