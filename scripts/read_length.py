#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute statistics over read lengths
"""


import gzip
import numpy

lengths = []

for fastq_path in snakemake.input:
    print(fastq_path)
    with gzip.open(fastq_path, 'rb') as fastq:
        for line in fastq:
            if not line.startswith((b"@", b"+")):
                lengths.append(len(line[:-1]))

sample = snakemake.wildcards["sample"]
lengths = numpy.array(lengths)
mean = lengths.mean()
std = lengths.std()

with open(snakemake.output[0], "w") as outlength:
    outlength.write(f"{sample}\t{mean}\t{std}\n")
