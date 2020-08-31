#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This script parses and formats picard collectalignmentsummarymetrics
information
"""

import pandas

data = pandas.read_csv(
    snakemake.input[0],
    sep="\t",
    header=0,
    index_col=None,
    skip_blank_lines=True,
    comment="#",
    dtype=str
)
data = data.loc[0]

headers = "\t".join(
    ["Sample_id", "Mean_insert_size", "Std_insert_size", "Orientation"]
)
results = "\t".join([
    snakemake.wildcards["sample"],
    data["MEAN_INSERT_SIZE"],
    data["STANDARD_DEVIATION"],
    data["PAIR_ORIENTATION"]
])
with open(snakemake.output[0], "w") as outfile:
    outfile.write(f"{results}\n")
