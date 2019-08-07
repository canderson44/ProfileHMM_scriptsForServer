#!/usr/bin/env python

'''
for each cell line (2_B01, 3_C01, 4_D01, splits filtered CCS sequences into
training and testing sets.
Creates a training and a testing file, fasta formatted.
Uses 10000 sequences for testing, 10000 for training. All randomly selected.
'''
import numpy as np
import random as rd

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/'
cells = ['2_B01', '3_C01', '4_D01']
num_train = 1000
num_test = 1000
for i in np.arange(3):
    cell = cells[i]
    input_filename = path_stub + cell + '_filtered_ccs.fasta'
    training_filename = path_stub + 'HMMConverter/cell_lines/' + cell + '_trainingSequences.txt'
    testing_filename = path_stub + 'HMMConverter/cell_lines/' + cell + '_testingSequences.txt'
    all_ccs = [] #all filtered ccs. format [ [ZMW,ccs],...]
    with open(input_filename,'r') as filtered:
        lines = filtered.readlines()
        for line in lines:
            line = line.strip()
            line_list = line.split(',') #[ZMW number, ccs]
            all_ccs.append(line_list)
    ccs_to_use = rd.sample(all_ccs,k=(num_train + num_test)) #selection without replacement
    training_ccs = ccs_to_use[:num_train]
    test_ccs = ccs_to_use[num_train:]

    #now write test and train files, fasta format

    #training
    with open(training_filename, 'w') as train:
        line_index = 0
        len_training = len(training_ccs)
        for entry in training_ccs:
            zmw = entry[0]
            ccs = entry[1]
            name_string = '>' + zmw + '\n'
            ccs_string = ccs
            if line_index < len_training -1: #means we're not at end of list of CCSes
                ccs_string += '\n'
            train.write(name_string)
            train.write(ccs_string)
            line_index +=1

    #testing
    with open(testing_filename, 'w') as test:
        line_index = 0
        len_test = len(test_ccs)
        for entry in test_ccs:
            zmw = entry[0]
            ccs = entry[1]
            name_string = '>' + zmw + '\n'
            ccs_string = ccs
            if line_index < len_test -1: #means we're not at end of list of CCSes
                ccs_string += '\n'
            test.write(name_string)
            test.write(ccs_string)
            line_index +=1