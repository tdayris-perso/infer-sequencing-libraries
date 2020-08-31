#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This script merges result metrics based on sample_id
"""

import pandas  # Handle tables and merge
import logging  # Traces and loggings

logging.basicConfig(
        filename=snakemake.log[0], filemode="w", level=logging.DEBUG
    )

data = None
for table in snakemake.input:
    frame = pandas.read_csv(table, sep="\t", header=0)
    frame.set_index("Sample_id", inplace=True)
    logging.debug(frame.head())
    try:
        data = pandas.merge(
            data,
            frame,
            left_index=True,
            right_index=True,
            how="outer",
            validate="1:1"
        )
    except TypeError:
        data = frame

data.to_csv(snakemake.output[0], sep="\t", index=True, header=True)
