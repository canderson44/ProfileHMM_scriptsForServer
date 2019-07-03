#!/usr/bin/env python

#created on: 2019-07-02
#author    : Catherine Anderson
#contact   : canderson44@wisc.edu


#uses generate_emission_probs.py and generate_transition_probs.py
#locations of above files:
    #/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/python_virtual_envs/envir0/lib/python3.6/site-packages

#for each cell type (2_B01, 3_D01, 4_D01), creates file called {cell}_hmm_initial.txt
#for ToPS program.
#The file specifies the initial parameters of the HMM model for use in the Baum-Welch
#parameter estimation process (done through ToPS)
#location of files written:

import generate_emission_probs as ep
import generate_transition_probs as tp
import numpy as np

transition_probs_list = tp.get_transition_probs()
emission_probs_list = ep.get_all_cells_emissions()

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/ToPS_files/'

cells = ['2_B01', '3_C01', '4_D01']

for i in np.arange(3):
    cell = cells[i]
    write_filename = path_stub + cell + '_hmm_initial.txt'

    #state names list

    match_names = ["M"+str(n) for n in np.arange(141)]
    match_rev_names = ["Mr"+str(n) for n in np.arange(141)]
    insert_names = ["I"+str(n) for n in np.arange(141)]
    insert_rev_names = ["I"+str(n) for n in np.arange(141)]
    delete_names = ["D"+str(n) for n in np.arange(139)]
    delete_rev_names = ["Dr"+str(n) for n in np.arange(139)]
    #start, end, junk states, RNAinsert states
    misc_names = ['START','END', 'IS', 'ISr', 'DS', 'DSr', 'R', 'Rr', 'IR',
                  'IRr', 'DR', 'DRr']
    #combine into one list
    state_names_list = match_names + match_rev_names + insert_names + insert_rev_names + delete_names + delete_rev_names + misc_names
    states_string = ""
    for entry in state_names_list:
        if states_string != "":
            states_string = states_string + ', ' + "\"" + entry + "\""
        else: #first entry
            states_string = "\"" + entry + "\""
 #   print(states_string)

    #format transitions: store in list
    #key is tuple:
        #P(first|second)
    transition_strings = []
    num_items = len(transition_probs_list[i].values())
    for key,value in transition_probs_list[i].items():
        first_part = "\"" + key[0] + "\""
        prob_string = first_part.ljust(10) + "| \"" + key[1] +"\": " + str(value)
        if len(transition_strings) < num_items - 1 : #not adding last element
            prob_string = prob_string + ";"
        transition_strings.append(prob_string)


    #format emissions: store in list
    #key is state, value is probability dictionary for emissions
    #so format will be:
        #   "A"    | "<key>" : <value['A']>;
        # and so on for C, G, T
    nucleotides = ["A", "C", "G", "T"]
    emission_strings = []
    num_items = len(emission_probs_list[i].values())
    for key,value in emission_probs_list[i].items():
        current_string = ""
        for nuc in nucleotides:
            first_part = "\"" + nuc + "\""
            current_string = first_part.ljust(5) + "| \"" + key +"\": " + str(value[nuc])
            if len(emission_strings) < num_items*4 - 1: #not adding last element:
                current_string = current_string + ";"
            emission_strings.append(current_string)

 #   print()


    #write file
    with open(write_filename, 'w') as f:
        f.write("# initial parameters for cell: " + cell + '\n')
        f.write("model_name = \"HiddenMarkovModel\" \n")
        f.write("state_names= (" + states_string + ') \n')
        f.write("observation_symbols= (\"A\", \"C\", \"G\", \"T\" ) \n")

        #transitions
        f.write("# transition probabilities \n")
        f.write("transititions = (" + transition_strings[0] + "\n")
        for n in np.arange(1,len(transition_strings)-1): #all but first and last entry
            f.write("".rjust(17) + transition_strings[n] + '\n')
        f.write("".rjust(17) + transition_strings[-1] + " )\n")

        #emissions
        f.write("# emission probabilities \n")
        f.write("emission_probabilities = (" + emission_strings[0] +"\n")
        for n in np.arange(1, len(emission_strings) - 1):  # all but first and last
            f.write("".rjust(26) + emission_strings[n] + '\n')
        f.write("".rjust(26) + emission_strings[-1] + ")\n")

        #initial probs
        f.write("initial_probabilities= (\"START\": 1.0)")