#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 19:18:03 2017

@author: michal kusibab, contribution paulina dziubanska-kusibab
"""

from Bio import SeqIO

genome_path='/path/to/human/genome/in/format.fasta'
bed_file_path='/path/to/BLISS/data/in/format.bed'
output_path='/path/to/output/file/file.fasta'

hg_dict = None
hg_dict = SeqIO.to_dict(SeqIO.parse(open(genome_path), 'fasta'))

left_shift=2
right_shift=1

chr=[str(x) for x in range(1,23)]+['X', 'Y']

output = open(output_path, 'w')
for line in open(bed_file_path):
   [id, start, end, cnt]= [line.split()[i] for i in [0,1,2,5]]
   if id in hg_dict and id in chr:
   #if id in hg_dict:
      output.write('>' + id + '_start:' + start + '_end:' + end + '_cnt:' + cnt +'\n' 
                   + str(hg_dict[id].seq[int(start)-left_shift:int(end)+1+right_shift]) + '\n')
   else:
      print id + ' was not found'
output.close()
