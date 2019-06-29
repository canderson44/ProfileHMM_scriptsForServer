#!/usr/bin/env python                                                           

#NOT CURRENT: for each cell line, takes first 1000 lines of filtered CCS csv 
#                           (found in {cell}_firstThousandLines.csv) 
#CURRENT: using all filtered CCSs becasue with 1000 ccs only got a couple successful strings written
#outputs file with each CCS annotated with the 5' barcode and appropriate 
#     adapter and 3' barcode
#annotation only added to file if adapter present in the annotation         
#uses Annotate_CCS.py in envir0 virtural environment
#between a region and its reverse complement, only includes one with higher score
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
    #    input_filename = path_stub + cell + '_firstThousandLines.csv'
    input_filename = path_stub + cell + '_filtered_ccs.csv'
    output_filename = path_stub + 'annotated_ccs/' + cell + '_annotatedCCS_findAdapter.txt'
    to_write_list = []
    with open(input_filename) as input:
        for line in input:
            if ',' in line:
                strip_line = line.strip()
                split_list = strip_line.split(',')
                to_write_list.append(split_list[1])                 
    with open(output_filename, 'w') as output:
        for line in to_write_list:
            ccs = line
            annotated_ccs = str(ac.annotate_seq(ccs, ref_list))
            if 'A' in annotated_ccs:
                output.write(annotated_ccs + '\n')
                output.write('\n')
