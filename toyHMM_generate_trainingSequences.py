#!/usr/bin/env python

'''
Writes sample sequences to file in desired ToPS format:
    "line starts with the sequence name, followed by a
    colon and the input sequence (in this section we
    assume the second option)"
                        --ToPS tutorial
Sequences follow emission probs
'''

#first: toy HMM Simple

''' 
TRANSITIONS
start -> M0: 1
M0 -> M1: 0.95
M0 -> D0: 0.05
D0 -> M2: 1
M1 -> M2: 1
M2 -> end: 1

EMISSIONS
M0: 0.95 A, rest 0.05
M1: 0.95 C, "   "
M2: 0.95 G, "   "
'''



import numpy as np
import random as rd

num_sequences = 1000
sequence_names = [str(n) for n in np.arange(num_sequences)]
#put into sample space 95 of .95 prob emission. rest, put in five.
ninetyfive_As = ['A'] * 95
ninetyfive_Cs = ['C'] * 95
ninetyfive_Gs = ['G'] * 95
five_As= ['A'] * 5
five_Cs = ['C'] * 5
five_Gs = ['G'] * 5
five_Ts = ['T'] * 5

m0_sample_space = ninetyfive_As + five_Cs + five_Gs + five_Ts
m1_sample_space = ninetyfive_Cs + five_As + five_Gs + five_Ts
m2_sample_space = ninetyfive_Gs + five_As + five_Cs + five_Ts

#transition sample spaces (only need m0's since others are all prob 1 to something

sample_space_transition_from_m0 = ['M1'] * 95
sample_space_transition_from_m0 += ['M2'] * 5 #represents going to D0, then M2

#now let's try make a list of sequences
all_sequences = []
all_transitions = []
for i in np.arange(num_sequences):
    #always start with M0
    char_list = []
    transition_list = ['M0']
    #first char
    char_list.append(rd.choice(m0_sample_space))
    #now transition
    transition = rd.choice(sample_space_transition_from_m0)
    transition_list.append(transition)

    while(transition != 'end'):
        if transition == 'M1':
            char_list.append(rd.choice(m1_sample_space))
            transition = 'M2'
        else: #transition to M2
            char_list.append(rd.choice(m2_sample_space))
            transition = 'end'
        transition_list.append(transition)
#    print("".join(char_list))
#    print(transition_list)

    all_sequences.append(" ".join(char_list))
    all_transitions.append("".join(transition))

output_name = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/HMMConverter/toyHMM_simple_sequences.fasta'
with open(output_name, 'w') as output:
    for index in np.arange(len(all_sequences)):
        label = '>' + sequence_names[index] + '\n'
        output.write(label)

        # don't make new line if this is the last sequence
        if index == len(all_sequences)-1:
            to_write_str = all_sequences[index]
        else:
            to_write_str = all_sequences[index] + '\n'
        output.write(to_write_str)


