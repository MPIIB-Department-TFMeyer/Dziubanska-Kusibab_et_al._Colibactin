#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 19:18:03 2017

@author: michal, 
edit: paulina (16.12.2019): 
# definition of output has changed, str(hg_dict[id].seq[int(start)-left_shift:int(end)+1+right_shift]) into: str(hg_dict[id].seq[int(start)-left_shift:int(end)+right_shift])
# cnt into strand;
# note: now for symetrical extraction left_shift = right_shift)
"""

from Bio import SeqIO

genome_path='/data_genome2/References/GenomeIndices/Human/human_g1k_v37_decoy_guideseq/human_g1k_v37_decoy.fasta'
# adapt this for each input file
bed_file_path='<BED input file containing DSB break points>'
output_path='<output FASTA file name>'

hg_dict = None
hg_dict = SeqIO.to_dict(SeqIO.parse(open(genome_path), 'fasta'))

left_shift=6
right_shift=6

chr=[str(x)for x in range(1,23)]+['X', 'Y']

output = open(output_path, 'w')
for line in open(bed_file_path):
   [id, start, end, strand]= [line.split()[i] for i in [0,1,2,5]]
   if id in hg_dict and id in chr:
   #if id in hg_dict:
      output.write('>' + id + '_start:' + start + '_end:' + end + '_strand:' + strand +'\n' 
                   + str(hg_dict[id].seq[int(start)-left_shift:int(end)+right_shift]) + '\n')
   else:
      print id + ' was not found'
output.close()
