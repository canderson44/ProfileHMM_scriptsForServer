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
    match_names = ["M"+str(n) for n in np.arange(141)]
    print(match_names)
    print()
    match_rev_names = []
    insert_names = []
    delete_names = []
    misc_names = []

    with open(write_filename, 'w') as f:
        f.write("#initial parameters for cell" + cell + '\n')
        f.write("model_name = HiddenMarkovModel")
