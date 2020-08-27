#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute statistics over read qualitites
"""


import gzip
import numpy

qualities = numpy.array
phreds = numpy.array

phred33 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopq"

with gzip.open(snakemake.input[0], "r") as fastq:
    line_iter = iter(fastq)
    line = next(line_iter, None)
    while line is not None:
        if line.startswith("@"):
            line = next(line_iter, None)
            line = next(line_iter, None)
        phred = (33 if all(i in phred33 for i in line[:-1]) else 64)
        qualities.append(sum(ord(i) - phred for i in line[:-1]))
        phreds.append(phred)

mean = qualities.mean()
std = qualities.std()
phred = phreads.mean()

with open(snakemake.output[0], 'w') as outfile:
    outfile.write(f"{mean}\t{std}\t{phred}\n")
