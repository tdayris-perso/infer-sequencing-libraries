#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute statistics over read lengths
"""


import gzip
import numpy


with gzip.open(snakemake.input[0], 'r') as fastq:
    length = numpy.array(
        len(line[:-1])
        for line in fastq
        if not line.startswith(("@", "+"))
    )

sample = snakemake.wildcards["sample"]
mean = length.mean()
std = length.std()

with open(snakemake.output[0], "w") as outlength:
    outlength.wirte(f"{sample}\t{mean}\t{std}\n")
