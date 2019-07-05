#!/usr/bin/env python

'''
for each cell line (2_B01, 3_C01, 4_D01, splits filtered CCS sequences into
training and testing sets.
Creates a training and a testing file, formatted as required by ToPS;
"each line starts with the sequence name, followed by a colon and the
input sequence" (from ToPS tutorial's first page)
'''
import numpy as np

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/'
cells = ['2_B01', '3_C01', '4_D01']

for i in np.arange(3):
    cell = cells[i]
    input_filename = path_stub + cell + '_filtered_ccs.csv'
    training_filename = path_stub + 'ToPS_files/' + cell + '_training_tops.sequences'
    testing_filename = path_stub + 'ToPS_files/' + cell + '_testing_tops.sequences'
