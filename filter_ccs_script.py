#!/usr/bin/env python
#note: for this script should use my venv called envir0
import random as rd
import numpy as np
import CCS_Filter as cf

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/'
start_ccs_path_stub = path_stub + 'demultiplexed_full_bams/no_passes/'
cells = ['2_B01']#, '3_C01', '4_D01']
#TODO RESTORE CELLS TO ALL CELLS
destination_stub = path_stub + 'profile_hmm/'

#3p barcodes
barcodes = [cf.barcode_2, cf.barcode_3, cf.barcode_4]

for i in np.arange(len(barcodes)):
    cell = cells[i]
    start_ccs_path = start_ccs_path_stub + cell + '/extracted_ccs.csv'
    barcode = barcodes[i]
    filtered_filename_path = destination_stub + cell + '_filtered_ccs.csv'
    # print("start_ccs_path: " + start_ccs_path)
    # print("filtered_filename_path: " + filtered_filename_path)
    cf.filter_ccs(start_ccs_path, filtered_filename_path, barcode)
    
