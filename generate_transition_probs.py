#!/usr/bin/env python
#generates transition probs to and from insert, match (sequence), and RNA_insert positions
#Reverse and forward sequence directions accounted for.
#assume forward sequence is of format 5'barcode -> RNAinsert -> 3'barcode
#FORMAT:  dictionary of tuples: {(A,B):0.05, (C,E):0.5}
#in this example, we have prob of A given B is 0.05, prob of C given E is 0.5
import numpy as np

avg_lengths = [1641.4098732872096, 1539.26103329929, 1652.740349292529]
#subtract off 96 basepairs (5pbarcode + 3pbarcode)
avg_lens_adjusted = [n-96 for n in avg_lengths]
cells = ['2_B01', '3_C01', '4_D01']

#only difference in cell lines is prob of self-cycle-of  and from  RNA_insert


#generates transition probabilities for each cell line in a dictionary, then returns
#a list of those dictionaries.
# entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01

transition_dict_list = []  # entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01
for i in np.arange(len(cells)):
    adjusted_length = avg_lens_adjusted[i]
    cell = cells[i]
    new_dict = {}
    # transition from start
    # 50% chance forward sequence, 50% chance backward sequence
    new_dict[("M0", "START")] = 0.45 #90% of 0.5
    # insert state before forward sequence; for junk
    new_dict[("IS","START")] = 0.025 #5% OF 0.5
    #delete states; Start to M1, M2, ... , M24, RNA
    #total probability allotted to deletion transitions
    tot_startDel_prob = 0.05 * 0.5
    num_startDel_toStates = 25 #M1 through and including RNA
    for i in np.arange(1,26):
        if i == 25:
            to_str = "RNA"
        else:
            to_str = "M" + str(i)
        new_dict[(to_str, "START")] = tot_startDel_prob/num_startDel_toStates

    #RC from start transitions
    new_dict[("Mr0", "START")] = 0.45 #90% of 0.5
    # insert state before reverse sequence; for junk
    new_dict[("ISr","START")] = 0.025
    #delete transitions: Start to Mr1, Mr2, ..., Mr24, RNAr
    # total probability allotted to deletion transitions
    tot_startDel_prob = 0.05 * 0.5
    num_startDel_toStates = 25  # M1 through and including RNA
    for i in np.arange(1, 26):
        if i == 25:
            to_str = "RNAr"
        else:
            to_str = "Mr" + str(i)
        new_dict[(to_str, "START")] = tot_startDel_prob / num_startDel_toStates

    #transitions from starting junk states
    #forward
    new_dict[("IS", "IS")] = 0.5
    new_dict[("M0", "IS")] = 0.5
    #backward
    new_dict[("ISr", "ISr")] = 0.5
    new_dict[("Mr0", "ISr")] = 0.5

    #############
    #############
    #############
    ## FORWARD ##
    ## SEQUENCE ##
    #############
    #############
    #############
    # match before RNA insert: 5' barcode
    for m_index in np.arange(24):  # M0 to M23
        match_str = "M" + str(m_index)
        next_match = "M" + str(m_index +1)
        insert_str = "I" + str(m_index)
        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9
        #delete states: i+2 to RNA inclusive
        tot_del_prob = 0.05
        tot_delete_toStates = 25-(m_index+1)
        if tot_delete_toStates >0: #else at state 24; shouldn't happen
            for successor in np.arange(m_index+2, 25): #don't want immediate next, but two ahead and onward
                to_delete_str = "M" + str(successor)
                new_dict[(to_delete_str,match_str)] = tot_del_prob/tot_delete_toStates
            #also delete transition to RNA
            new_dict[("RNA",match_str)] = tot_del_prob/tot_delete_toStates


        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5
    # END MATCH before RNA-INSERT

    # M24: just before RNAinsert
    # no delete, no (Match, Match)
    # transitions: (R, M24), (I24, M24)
    new_dict[("I24", "M24")] = 0.05
    new_dict[("RNA", "M24")] = 0.95
    # I24 transitions
    new_dict[("I24", "I24")] = 0.5
    new_dict[("RNA", "I24")] = 0.5

    # transitions from RNAinsert
    selfCycle_prob = 1 - (1.0 / adjusted_length)
    new_dict[("RNA", "RNA")] = selfCycle_prob
    remaining_prob = 1 - selfCycle_prob
    five_percent = 0.05 * remaining_prob
    new_dict[("M25", "RNA")] = remaining_prob - five_percent - five_percent
    new_dict[("IRNA", "RNA")] = five_percent
    #delete transitions: M26, M27, ... , M95, End
    num_del_toStates = 71 #M26 to M95 and End
    tot_del_prob = five_percent
    for successor in np.arange(26,97):
        if successor == 96:
            to_delete_str = "END"
        else:
            to_delete_str = "M" + str(successor)
        new_dict[(to_delete_str, "RNA")] = tot_del_prob/num_del_toStates

    # transitions from IR
    new_dict[("IRNA", "IRNA")] = 0.5
    new_dict[("M25", "IRNA")] = 0.5

    # after RNA: M25 to M95 and END: 3' barcode
    for m_index in np.arange(25, 95):  # excludes last match state, M95
        match_str = "M" + str(m_index)
        next_match = "M" + str(m_index + 1)
        insert_str = "I" + str(m_index)

        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9
        #Deletions: transitions to all succeeding states starting with M_index + 2, to END
        tot_del_prob = 0.05
        num_del_toStates = 96 - (m_index + 1)
        if num_del_toStates >0:
            for successor in np.arange(m_index+2, 96): #transitions up to including M95
                to_delete_str = "M" + str(successor)
                new_dict[(to_delete_str, match_str)] = tot_del_prob/num_del_toStates
            #also transition to end state
            new_dict[("END", match_str)] = tot_del_prob/num_del_toStates

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

    # transitions from M95 (last match state)
    # no delete, no (Match, Match)
    # transitions: (END, M95), (I95, M95)
    # transition from match state M
    new_dict[("I95", "M95")] = 0.05
    new_dict[("END", "M95")] = 0.95
    # I transitions
    new_dict[("I95", "I95")] = 0.5
    new_dict[("END", "I95")] = 0.5

    ##############
    ##############
    ##############
    ## REVERSE ###
    # COMPLEMENT #
    ##############
    ##############
    ##############

    # match before RNA insert: 3' barcode RC
    for m_index in np.arange(70):  # M0 to Mr69; excludes last state M70
        match_str = "Mr" + str(m_index)
        next_match = "Mr" + str(m_index + 1)
        insert_str = "Ir" + str(m_index)
        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9
        #deletions: transition to M[m_index+2] up to RNAr inclusive
        tot_del_prob = 0.05
        tot_delete_toStates = 71 -(m_index + 1)
        if tot_delete_toStates >0:
            for successor in np.arange(m_index+2,71): #transitions up to including Mr70
                to_delete_str = "Mr" + str(successor)
                new_dict[(to_delete_str, match_str)] = tot_del_prob/tot_delete_toStates
            #also deletion transition to RNAr
            new_dict[("RNAr", match_str)] = tot_del_prob/tot_delete_toStates

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5
    # END MATCH before RNA-INSERT

    # Mr70: just before RNAinsert
    # no delete, no (Match, Match)
    # transitions: (Rr, Mr70), (Ir95, Mr70)
    # transition from match state M
    new_dict[("Ir70", "Mr70")] = 0.05
    new_dict[("RNAr", "Mr70")] = 0.95
    # I24 transitions
    new_dict[("Ir70", "Ir70")] = 0.5
    new_dict[("RNAr", "Ir70")] = 0.5

    # Transitions from RNAr
    new_dict[("RNAr", "RNAr")] = selfCycle_prob
    remaining_prob = 1 - selfCycle_prob
    five_percent = 0.05 * remaining_prob
    tot_del_prob = five_percent
    #next match: M71
    new_dict[("Mr71", "RNAr")] = remaining_prob - five_percent - five_percent
    new_dict[("IRNAr", "RNAr")] = five_percent
    #RNAr delete transitions: to Mr72, Mr73, ... , Mr95, END len 25
    tot_delete_toStates = 25
    for successor in np.arange(72,97):
        if successor == 96:
            to_delete_str = "END"
        else:
            to_delete_str = "Mr" + str(successor)
        new_dict[(to_delete_str, "RNAr")] = tot_del_prob/tot_delete_toStates

    # transitions from IRNAr
    new_dict[("IRNAr", "IRNAr")] = 0.5
    new_dict[("Mr71", "IRNAr")] = 0.5

    # after RNAr
    for m_index in np.arange(71, 95):  # excludes last match state Mr95
        match_str = "Mr" + str(m_index)
        next_match = "Mr" + str(m_index + 1)
        insert_str = "Ir" + str(m_index)

        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9
        tot_del_prob = 0.05
        num_del_toStates = 96 - (m_index + 1)
        #deletions: to M[m_index + 2] up to including END
        if num_del_toStates >0:
            for successor in np.arange(m_index+2, 96): #Mr[m_index + 2] to Mr95
                to_delete_str = "Mr" + str(successor)
                new_dict[(to_delete_str,match_str)] = tot_del_prob/num_del_toStates
            #also deletion transition to END
            new_dict[("END", match_str)] = tot_del_prob/num_del_toStates

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

    # transitions from Mr95. Last match state
    # no delete, no (Match, Match)
    # transitions: (END, M95), (I95, M95)
    new_dict[("Ir95", "Mr95")] = 0.05
    new_dict[("END", "Mr95")] = 0.95
    # I140 transitions
    new_dict[("Ir95", "Ir95")] = 0.5
    new_dict[("END", "Ir95")] = 0.5
    transition_dict_list.append(new_dict)

#getter
#returns list of transition dictionaries, one per cell
## entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01
def get_transition_probs():
    return transition_dict_list

#getter
#returns 2B01's transition probability dictionary
def get_2B01_transitions():
    return transition_dict_list[0]

#getter
#returns 3C01's transition probability dictionary
def get_3C01_transitions():
    return transition_dict_list[1]

#getter
#returns 4D01's transition probability dictionary
def get_4D01_transitions():
    return transition_dict_list[2]


#for testing purposes:
# for key,value in get_2B01_transitions().items():
#     print(str(key) + ':' + str(value))