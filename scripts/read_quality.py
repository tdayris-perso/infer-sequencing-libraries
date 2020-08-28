#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute statistics over read qualitites
"""


import gzip
import numpy

qualities = []


for fastq_path in snakemake.input:
    with gzip.open(fastq_path, "rb") as fastq:
        print(fastq_path)
        line_iter = iter(fastq)
        line = next(line_iter, None)
        while line is not None:
            if line.startswith(b"@"):
                line = next(line_iter, None)
                line = next(line_iter, None)
                line = next(line_iter, None)
            qualities.append(sum(ord(i) - 33 for i in str(line[:-1])) / len(str(line[:-1])))
            line = next(line_iter, None)
            line = next(line_iter, None)

sample = snakemake.wildcards["sample"]
qualities = numpy.array(qualities)
mean = qualities.mean()
std = qualities.std()

with open(snakemake.output[0], 'w') as outfile:
    outfile.write(f"{sample}\t{mean}\t{std}\n")
