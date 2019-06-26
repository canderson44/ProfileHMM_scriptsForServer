#!/usr/bin/env python

#for each cell line, writes first 10 lines of CCS csv file to a new file
#for now: input file is no_passes's extracted_ccs.csv
#TODO: once filtering done, use profile_hmm's filtered_ccs csv files

import numpy as np
path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/'
cells = ['2_B01', '3_C01', '4_D01']
for cell in cells:
    input_filename = path_stub + cell + '_filtered_ccs.csv'
    #3_C01_first10lines.csv
    output_filename = path_stub + cell + '_first10lines.csv'
    length_list = []
    with open(input_filename) as input:
        with open(output_filename, 'w+') as output:
            for line in input:
                write_line = line.strip()
                output.write(write_line + '\n')
                length_list.append(len(write_line.split(',')[1]))
            avg_length = sum(length_list)/len(length_list)
            output.write('avg length is:' + str(avg_length)) 
