#!/usr/bin/env python

#created on: 2019-07-02
#author    : Catherine Anderson
#contact   : canderson44@wisc.edu

'''
uses generate_emission_probs.py and generate_transition_probs.py
locations of above files:
    /tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/demultiplexed_full_bams/bash_scripts/

for each cell type (2_B01, 3_D01, 4_D01), creates file called {cell}_hmm_initial.xml
for ToPS program.
The file specifies the initial parameters of the HMM model for use in the Baum-Welch
parameter estimation process (done through HMMCompiler)
location of files written:
    /tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/HMMConverter/cell_lines
'''

import generate_emission_probs as ep
import generate_transition_probs as tp
import numpy as np

transition_probs_list = tp.get_transition_probs()
emission_probs_list = ep.get_all_cells_emissions()

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/HMMConverter/cell_lines/'

cells = ['2_B01', '3_C01', '4_D01']

INDENT = '        '#eight spaces
FOUR_SPACES = '    '#just four spaces
nucleotides = ['A','C','G','T']
for i in np.arange(len(cells)):
    cell = cells[i]
    write_filename = path_stub + cell + '_hmm_initial.xml'
    emission_filename = path_stub + cell + '_emissions_initial.txt'
    seq_filename = cell + '_trainingSequences.txt'

    #state names list
    match_names = ["M"+str(n) for n in np.arange(96)]
    match_rev_names = ["Mr"+str(n) for n in np.arange(96)]
    insert_names = ["I"+str(n) for n in np.arange(96)]
    insert_rev_names = ["Ir"+str(n) for n in np.arange(96)]
    #start, end, junk states, RNAinsert states
    misc_names = ['IS', 'ISr', 'RNA', 'RNAr', 'IRNA',
                  'IRNAr']
    #combine into one list. Want start first, end last
    state_names_list = ["START"] + match_names + match_rev_names + insert_names + insert_rev_names + misc_names + ["END"]
    emitting_states = match_names + match_rev_names + insert_names + insert_rev_names + misc_names
    # first: get transitions
    # format transitions: store in list
    # key is tuple:
    # P(first|second)
    # value is prob
    transition_dict = dict()
    # format of transition_dict:
    # key is from name, value is list of (to_name, probability) tuples
    num_items = len(transition_probs_list[i].values())
    for key, value in transition_probs_list[i].items():
        from_state_name = key[1]
        to_state_name = key[0]
        if from_state_name in transition_dict:  # already a list for this state
            transition_dict[from_state_name].append((to_state_name, value))
        else:  # no list yet
            transition_dict[from_state_name] = [(to_state_name, value)]


    #get emissions
    emissions_dict = emission_probs_list[i]
    #format for emissions_dict:
    # key is state name, value is probability dictionary for emissions
    ##################
    ##################
    ##################
    ##################
    ##################
    # Emission file  #
    ##################
    ##################
    ##################
    ##################
    ##################
    ##################
    ##################
    with open(emission_filename, 'w') as output:
        for index in np.arange(len(emitting_states)):
            name = emitting_states[index]
            id = "EP." + str(index)
            first_line = id + FOUR_SPACES + name + FOUR_SPACES + "1\n"
            output.write(first_line)
            emission_lines = ''
            for nuc in nucleotides:
                prob = emissions_dict[name][nuc]
                emission_lines += nuc + ' ' + str(prob) + '\n'
            if index<len(emitting_states)-1: #not last one
                emission_lines += '\n'
            output.write(emission_lines)








    ##################
    ##################
    ##################
    ##################
    ##################
    # Main input file#
    ##################
    ##################
    ##################
    ##################
    ##################
    ##################
    ##################



    #start writing file
    with open(write_filename,'w') as output:
        first_lines = '<?xml version="1.0"?>\n<HMMConverter>\n\n<model>\n'
        model_type_tag = INDENT + '<Model_Type name=\"' + cell + 'HMM\"/>\n\n'
        alphabet_tag = INDENT + '<Alphabets set=\"ACGT\" />\n\n'
        emissions_tag = INDENT + '<Emission_Probs id=\"EP\" size=\"390\" file=\"' + emission_filename + '\"/>\n\n'
        states_beg_tag = INDENT + '<States>\n'
        output.write(first_lines)
        output.write(model_type_tag)
        output.write(alphabet_tag)
        output.write(emissions_tag)
        output.write(states_beg_tag)


        ########
        #States#
        ########

        id_name_dict = dict() #key is name, value is id
        for index in np.arange(len(state_names_list)):
            #now list all the state tags
            this_indent = INDENT + INDENT
            name=state_names_list[index]
            id = "S." + str(index)
            #add id-name combo to dict
            id_name_dict[name] = id
            if name == 'START' or name=='END': #don't emmit any characters
                xdim = str(0)
            else:
                xdim = str(1)
            this_str = this_indent + "<State id=\"" + id + "\" name=\"" + name + "\" xdim=\"" + xdim + "\"/>\n"
            output.write(this_str)




        end_states_tag =INDENT+"</States>\n\n"
        start_transitions_tag = INDENT + "<Transitions>\n"
        output.write(end_states_tag)
        output.write(start_transitions_tag)

        #############
        #Transitions#
        #############
        #now write all transition tags
        for index in np.arange(len(state_names_list)): #iterate through each state
            from_indent = INDENT + INDENT
            to_indent = from_indent + INDENT
            from_name = state_names_list[index]
            from_id = id_name_dict[from_name]
            if from_name != "END":
                from_start_tag = from_indent + "<from idref=\"" + from_id + "\">\n"
                output.write(from_start_tag)
     #           to_tag_list = []
                for pair in transition_dict[from_name]: #format (to_name, prob)
                    to_name = pair[0]
                    to_id = id_name_dict[to_name]
                    prob = str(pair[1])
                    this_str = to_indent + "<to idref=\"" + to_id + "\" exp=\"" + prob + "\"/>\n"
                    output.write(this_str)
                end_from_tag = from_indent + "</from>\n"
                output.write(end_from_tag)
            else: #this is the END state
                from_tag = from_indent + "<from idref=\"" + from_id + "\"/>\n"
                output.write(from_tag)

        #end transition, end model
        end_transition_tag = INDENT + "</Transitions>\n"
        end_model_tag = "</model>\n\n"
        output.write(end_transition_tag)
        output.write(end_model_tag)

        ###################
        #Sequence Analysis#
        ###################

        start_seqAnalysis_tag = "<sequence_analysis>\n"
        output.write(start_seqAnalysis_tag)
        start_param_tag = INDENT + "<parameter_training>\n"
        output.write(start_param_tag)
        #max vol, max iter, threshold for Baum Welch totally arbitrary
        maxVol = 1000000
        maxIter = 100
        threshold = 0.00001
        input_tag = INDENT + INDENT + "<input_files SeqFile=\"" + seq_filename + "\"/>\n"
        alg_tag = INDENT + INDENT + "<algorithm alg=\"0\" MaxVolume=\"" + str(maxVol) + "\" Maxiter=\"" \
                  + str(maxIter) + "\" threshold=\"" + str(threshold) + '\"/>\n'
        output_tag = INDENT + INDENT + "<output_files XMLFile=\"" + cell + "_trainedHMM.xml\""
        output_tag += " EProbFile=\"" + cell + "_trainedEmissions.txt\"/>\n"
        output.write(input_tag)
        output.write(alg_tag)
        output.write(output_tag)

        #wrap up parameter training, sequence analysis, and HMMConverter tags
        end_param_tag = INDENT + "</parameter_training>\n"
        end_seqAnalysis_tag = "</sequence_analysis>\n\n"
        end_HMMConverter_tag = "</HMMConverter>\n"
        output.write(end_param_tag)
        output.write(end_seqAnalysis_tag)
        output.write(end_HMMConverter_tag)

    # if cell == '2_B01':
    #     print (id_name_dict)


