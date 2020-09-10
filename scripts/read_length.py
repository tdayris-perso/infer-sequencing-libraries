#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute statistics over read lengths
"""


import gzip
import numpy
import logging

logging.basicConfig(
    filename=snakemake.log[0], filemode="w", level=logging.DEBUG
)

lengths = []

for fastq_path in snakemake.input:
    logging.debug(f"Processing: {fastq_path}")
    with gzip.open(fastq_path, 'rb') as fastq:
        for line in fastq:
            if not line.startswith((b"@", b"+")):
                lengths.append(len(line[:-1]))

sample = snakemake.wildcards["sample"]
lengths = numpy.array(lengths)
mean = lengths.mean()
std = lengths.std()
results = f"{sample}\t{mean}\t{std}\n"
logging.debug(results)

with open(snakemake.output[0], "w") as outlength:
    outlength.write(results)
