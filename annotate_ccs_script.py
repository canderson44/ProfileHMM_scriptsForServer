#!/usr/bin/env python

#for each cell line, takes first 10 lines of filtered CCS csv 
#                           (found in {cell}_first10lines.csv)
#outputs file with each CCS annotated with the 5' barcode and appropriate 
#     adapter and 3' barcode 
#uses Annotate_CCS.py in envir0 virtural environment

import Annotate_CCS.py as ac
import numpy as np

cells = ['2_B01', '3_C01', '4_D01']
barcode_chars = ['2','3','4']
barcodes3p = [ac.barcode2, ac.barcode3, ac.barcode4]
adapters = [ac.two_adapter, ac.three_adapter, ac.four_adapter]
path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/'
for i in np.arange(3):
    cell = cells[i]
    ref_list = [[ac.fivePBarcode,'5'], [barcodes3p[i],barcode_chars[i]], [adapters[i],'A']]
    input_filename = path_stub + cell + '_first10lines.csv'
    output_filename = path_stub + 'annotated_ccs/' + cell + '_annotatedCCS.txt'
    with open(input_filename) as input:
        with open(output_filename, 'w+') as output:
            for line in input:
                line = line.strip()
                ccs = line.split(',')[1]
                output.write(ac.annotate_seq(ccs,ref_list) + '\n')
                output.write('\n')


            
