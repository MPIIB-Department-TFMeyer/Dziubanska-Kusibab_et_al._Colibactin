#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:18:19 2017

@author: hilmar
"""

import gzip

ifile_name = "gencode.v28lift37.transcripts.fa.gz"
ofile_name = "gencode.v28lift37.transcript.anno.txt"

with open(ofile_name, "w") as ofile:
    ofile.write("transcript\tgene\tOTTgene\tOTTtx\ttxSymbol\tGeneSymbol\tlength\tTxType\n")
    with gzip.open(ifile_name,"r") as ifile:
        for line in ifile.readlines():
            if line[0] != ">":
                continue
            fields = line[1:].strip("\n\r|").split("|")
            ofile.write("\t".join(fields) + "\n")
