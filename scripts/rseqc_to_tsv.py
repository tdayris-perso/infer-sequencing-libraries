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
    "Protocol"
])

with open(snakemake.input[0]) as rseqc:
    fractions = [
        line.split(":")[-1][:-1]
        for line in rseqc
        if not (line == "\n" or line.startswith("This is"))
    ]
logging.debug(f"{snakemake.wildcards.sample} processed")


orientation = None

if all(0.40 < frac < 0.60 for frac in map(float, fractions[-2:])):
    orientation = "Unstranded"
elif all(0.10 < frac < 0.90 for frac in map(float, fractions[-2:])):
    orientation = "Stranded"
else:
    orientation = "Unknown"

content = "\t".join([snakemake.wildcards.sample] + fractions + [orientation])

with open(snakemake.output[0], 'w') as outtsv:
    outtsv.write("\n".join([titles, content]))
    logging.debug("Process over")
