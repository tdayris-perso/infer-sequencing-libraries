#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This script parses and formats picard collectalignmentsummarymetrics
information
"""

import pandas
import logging

logging.basicConfig(
    filename=snakemake.log[0], filemode="w", level=logging.DEBUG
)

data = pandas.read_csv(
    snakemake.input[0],
    sep="\t",
    header=0,
    index_col=0,
    skip_blank_lines=True,
    comment="#"
)
logging.debug(data.head())


# Guessing library
library = "Unknown"
total_fr = data.loc["FIRST_OF_PAIR", "TOTAL_READS"]
total_sr = data.loc["SECOND_OF_PAIR", "TOTAL_READS"]

if (total_fr - total_fr * 0.1) < total_sr < (total_fr + total_fr * 0.1):
    library = "Unstranded"
elif (total_fr - total_fr * 0.9) < total_sr < (total_fr + total_fr * 0.9):
    library = "Stranded"


# Printing results
headers = "\t".join(
    ["Sample_id", "Upstream_reads", "Downstream_reads", "Library"]
)
results = "\t".join([
    snakemake.wildcards["sample"],
    str(total_fr),
    str(total_sr),
    library
])
logging.head(results)
with open(snakemake.output[0], "w") as outfile:
    outfile.write(f"{results}\n")
    logging.debug("Process over.")
