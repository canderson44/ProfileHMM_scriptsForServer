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
start -> I0: 0.05
start -> M0: 0.95
I0 -> I0: 0.5
I0 -> M0: 0.5
M0 -> M1: 0.95
M0 -> M2: 0.05
M1 -> M2: 1
M2 -> I2: 1
I2 -> I2: 0.80
I2 -> END: 0.20

EMISSIONS
M0: 0.95 A, rest 0.05
M1: 0.95 C, "   "
M2: 0.95 G, "   "
I0 and I2 0.25 for each
'''



import numpy as np
import random as rd

num_sequences = 1000
sequence_names = [str(n) for n in np.arange(num_sequences)]
#put into sample space 95 of .95 prob emission. rest, put in five.
ninetyfive_As = ['a'] * 95
ninetyfive_Cs = ['c'] * 95
ninetyfive_Gs = ['g'] * 95
nucleotides = ['A', 'C', 'G', 'T']
five_As= ['a'] * 5
five_Cs = ['c'] * 5
five_Gs = ['g'] * 5
five_Ts = ['t'] * 5

m0_sample_space = ninetyfive_As + five_Cs + five_Gs + five_Ts
m1_sample_space = ninetyfive_Cs + five_As + five_Gs + five_Ts
m2_sample_space = ninetyfive_Gs + five_As + five_Cs + five_Ts

#transition sample spaces (only need m0's since others are all prob 1 to something
sample_space_transition_from_I0 = ['I0', 'M0']
sample_space_transition_from_I2 = ['I2'] * 80
sample_space_transition_from_I2 += ['END'] * 20
sample_space_transition_from_m0 = ['M1'] * 95
sample_space_transition_from_m0 += ['M2'] * 5 #represents going to D0, then M2
sample_space_startState = ['M0'] * 95
sample_space_startState += ['I0'] * 5
#now let's try make a list of sequences
converter_all_sequences = []
tops_all_sequences = []
all_transitions = []
transition = rd.choice(sample_space_startState)
for i in np.arange(num_sequences):
    char_list = []
    transition_list = []
    while(transition != 'end'):
        if transition == 'I0':
            char_list.append(rd.choice(nucleotides))
            transition= rd.choice(sample_space_transition_from_I0)
        elif transition == 'M0':
            char_list.append(rd.choice(m0_sample_space))
            transition = rd.choice(sample_space_transition_from_m0)
        elif transition == 'M1':
            char_list.append(rd.choice(m1_sample_space))
            transition = 'M2'
        elif transition == 'M2': #transition to M2
            char_list.append(rd.choice(m2_sample_space))
            transition = 'I2'
        else: #transition is I2
            char_list.append(rd.choice(nucleotides))
            transition = rd.choice(sample_space_transition_from_I2)
        transition_list.append(transition)
#    print("".join(char_list))
#    print(transition_list)

    converter_all_sequences.append("".join(char_list)) #FOR HMMCONVERTER
    tops_all_sequences.append(" ".join(char_list))
    all_transitions.append("".join(transition))

converter_output_name = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/HMMConverter/toyHMM_simple_sequences.txt'
tops_output_name = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/ToPS/toyHMM_simple.sequences'
with open(tops_output_name, 'w') as tops:
    with open(converter_all_sequences,'w') as converter:
        for index in np.arange(len(tops_all_sequences)):
            converter_label = '>' + sequence_names[index] + '\n' #FOR HMMCONVERTER
            converter.write(converter_label)                         #FOR HMMCONVERTER
            tops_label= sequence_names[index] + ': '

            # don't make new line if this is the last sequence
            if index == len(tops_all_sequences)-1:
                converter_to_write_str = converter_all_sequences[index] #FOR HMMCONVERTER
                tops_to_write_str = tops_label + tops_all_sequences[index]
            else:
                converter_to_write_str = converter_all_sequences[index] + '\n' #FOR HMMCONVERTER
                tops_to_write_string = tops_label + tops_all_sequences[index] + '\n'
            tops.write(tops_to_write_str)
            converter.write(converter_to_write_str)


