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

cells = ['2_B01', '3_D01', '4_D01']

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
    print(states_string)

    #format transitions: store in list
    transition_strings = []
    for key,value in transition_probs_list[i]:
        print("key is: " + key)
        print("value is: " + value)

    #write file
    with open(write_filename, 'w') as f:
        f.write("# initial parameters for cell: " + cell + '\n')
        f.write("model_name = HiddenMarkovModel + \n")
        f.write("state_names= (" + states_string + ') \n')
        f.write("observation_symbols= (\"A\", \"C\", \"G\", \"T\" ) \n")
        f.write("# transition probabilitites \n")
        f.write("transititions = (" + ') \n')