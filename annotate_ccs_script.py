#!/usr/bin/env python

#for each cell line, takes first 10 lines of filtered CCS csv 
#                           (found in {cell}_first10lines.csv)
#outputs file with each CCS annotated with the 5' barcode and appropriate 
#     adapter and 3' barcode 
#uses Annotate_CCS.py in envir0 virtural environment

import Annotate_CCS as ac
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
    to_write_list = []
    with open(input_filename) as input:
        for line in input:
            if ',' in line:
                strip_line = line.strip()
                split_list = strip_line.split(',')
#                print(split_list)
                to_write_list.append(split_list[1])
#               print("successfully appended")
#            print()
    with open(output_filename, 'w') as output:
        for line in to_write_list:
            ccs = line
            annotated_ccs = str(ac.annotate_seq(ccs, ref_list))
#            print("annotated_ccs:",annotated_ccs)
#            print()
            #output.write("testing 1234     \n")
            output.write(annotated_ccs + '\n')


            
