#!/usr/bin/env python
#generates transition probs to and from delete, insert, match (sequence), and RNA_insert positions
#for now, just a high prob: 0.85
# TODO add reverse complement!!!!!!
#FORMAT:  dictionary of tuples: {(A,B):0.05, (C,E):0.5}
#in this example, we have prob of A given B is 0.05, prob of C given E is 0.5
import numpy as np

avg_lengths = [1641.4098732872096, 1539.26103329929, 1652.740349292529]
#subtract off 141 basepairs (5pbarcode + 3pbarcode, adapter)
avg_lens_adjusted = [n-141 for n in avg_lengths]
cells = ['2_B01', '3_C01', '4_D01']

#only difference in cell lines is prob of self-cycle-of  and from  RNA_insert


#generates transition probabilities for each cell line in a dictionary, then returns
#a list of those dictionaries.
# entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01

transition_dict_list = []  # entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01
match_b4_RNAinsert = 25  # just 5p barcode
match_after_RNAinsert = 116  # adapter + 3p barcode
for i in np.arange(3):
    cell = cells[i]
    new_dict = {}
    # transition from start
    # 50% chance forward sequence, 50% chance backward sequence
    new_dict[("M0", "START")] = 0.45 #90% of 0.5
    # insert state before forward sequence; for junk
    new_dict[("IS","START")] = 0.025 #5% OF 0.5
    #delete state; skip first element of 5p barcode
    new_dict[("DS","START")] = 0.025 #5% OF 0.5
    new_dict[("Mr0", "START")] = 0.45 #90% of 0.5
    # insert state before reverse sequence; for junk
    new_dict[("ISr","START")] = 0.025
    # delete state; skip first element of 5p barcode
    new_dict[("DSr", "START")] = 0.025  # 5% OF 0.5

    #transitions from starting junk states
    #forward
    new_dict[("IS", "IS")] = 0.5
    new_dict[("M0", "IS")] = 0.5
    #backward
    new_dict[("ISr", "ISr")] = 0.5
    new_dict[("Mr0", "ISr")] = 0.5

    #transitions from starting delete states
    #forward
    new_dict[("D0", "DS")] = 0.5
    new_dict[("M1", "DS")] = 0.5
    #backward
    new_dict[("Dr0", "DSr")] = 0.5
    new_dict[("Mr1", "DSr")] = 0.5
    #############
    #############
    #############
    ## FORWARD ##
    ## SEQUENCE ##
    #############
    #############
    #############
    # match before RNA insert
    for m_index in np.arange(24):  # M0 to M23
        match_str = "M" + str(m_index)
        next_match = "M" + str(m_index + 1)
        two_ahead_match = "M" + str(m_index + 2)
        del_str = "D" + str(m_index)
        next_del = "D" + str(m_index + 1)
        insert_str = "I" + str(m_index)
        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(del_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

        # transition from delete state
        if m_index < 23:
            new_dict[(next_del, del_str)] = 0.5
            new_dict[(two_ahead_match, del_str)] = 0.5
        else:  # no del->next_del transition
            # no del->match transition
            # only transition: d23 -> RNAinsert
            new_dict[("R", del_str)] = 1.0

    # END MATCH B4 RNAINSERT

    # M24: just before RNAinsert
    # no delete, no (Match, Match)
    # transitions: (R, M24), (I24, M24)
    # transition from match state M
    new_dict[("I24", "M24")] = 0.05
    new_dict[("R", "M24")] = 0.95
    # I24 transitions
    new_dict[("I24", "I24")] = 0.5
    new_dict[("R", "I24")] = 0.5

    # transitions from RNAinsert
    p = 1 - 1.0 / (avg_lens_adjusted[i])
    new_dict[("R", "R")] = p
    remaining_prob = 1 - p
    five_percent = 0.05 * remaining_prob
    new_dict[("M25", "R")] = remaining_prob - five_percent - five_percent
    new_dict[("DR", "R")] = five_percent
    new_dict[("IR", "R")] = five_percent
    # transitions from IR
    new_dict[("IR", "IR")] = 0.5
    new_dict[("M25", "IR")] = 0.5
    # transitions from DR
    new_dict[("D25", "DR")] = 0.5
    new_dict[("M26", "DR")] = 0.5

    # after RNA
    for m_index in np.arange(25, 140):  # excludes last match state, M140
        match_str = "M" + str(m_index)
        next_match = "M" + str(m_index + 1)
        two_ahead_match = "M" + str(m_index + 2)
        del_str = "D" + str(m_index)
        next_del = "D" + str(m_index + 1)
        insert_str = "I" + str(m_index)

        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(del_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

        # transition from delete state
        if m_index < 139:
            new_dict[(next_del, del_str)] = 0.5
            new_dict[(two_ahead_match, del_str)] = 0.5


        else:  # no del->next_del transition
            # no del->match transition
            # only transition: d139 -> END
            new_dict[("END", del_str)] = 1.0

    # transitions from M140. Last match state
    # no delete, no (Match, Match)
    # transitions: (END, M140), (I140, M140)
    # transition from match state M
    new_dict[("I140", "M140")] = 0.05
    new_dict[("END", "M140")] = 0.95
    # I24 transitions
    new_dict[("I140", "I140")] = 0.5
    new_dict[("END", "I140")] = 0.5

    ##############
    ##############
    ##############
    ## REVERSE ###
    # COMPLEMENT #
    ##############
    ##############
    ##############
    match_after_RNAinsert = 25  # just 5p barcode       Mr116 to Mr140
    match_b4_RNAinsert = 116  # adapter + 3p barcode   Mr0 to Mr115

    # match before RNA insert
    for m_index in np.arange(115):  # M0 to Mr114
        match_str = "Mr" + str(m_index)
        next_match = "Mr" + str(m_index + 1)
        two_ahead_match = "Mr" + str(m_index + 2)
        del_str = "Dr" + str(m_index)
        next_del = "Dr" + str(m_index + 1)
        insert_str = "Ir" + str(m_index)
        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(del_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

        # transition from delete state
        if m_index < 114:
            new_dict[(next_del, del_str)] = 0.5
            new_dict[(two_ahead_match, del_str)] = 0.5
        else:  # no del->next_del transition
            # no del->match transition
            # only transition: d114 -> RNAinsert
            new_dict[("Rr", del_str)] = 1.0

    # END MATCH B4 RNAINSERT

    # M115: just before RNAinsert
    # no delete, no (Match, Match)
    # transitions: (Rr, M115), (I114, M115)
    # transition from match state M
    new_dict[("Ir115", "Mr115")] = 0.05
    new_dict[("Rr", "Mr115")] = 0.95
    # I24 transitions
    new_dict[("Ir115", "Ir115")] = 0.5
    new_dict[("Rr", "Ir115")] = 0.5

    # Transitions from Rr
    new_dict[("Rr", "Rr")] = p
    remaining_prob = 1 - p
    five_percent = 0.05 * remaining_prob
    new_dict[("Mr116", "Rr")] = remaining_prob - five_percent - five_percent
    new_dict[("DRr", "Rr")] = five_percent
    new_dict[("IRr", "Rr")] = five_percent
    # transitions from IR
    new_dict[("IRr", "IRr")] = 0.5
    new_dict[("Mr116", "IRr")] = 0.5
    # transitions from DR
    new_dict[("Dr116", "DRr")] = 0.5
    new_dict[("Mr117", "DRr")] = 0.5

    # after RNA
    for m_index in np.arange(116, 140):  # excludes last match state, M140
        match_str = "Mr" + str(m_index)
        next_match = "Mr" + str(m_index + 1)
        two_ahead_match = "Mr" + str(m_index + 2)
        del_str = "Dr" + str(m_index)
        next_del = "Dr" + str(m_index + 1)
        insert_str = "Ir" + str(m_index)

        # transition from match state M
        new_dict[(insert_str, match_str)] = 0.05
        new_dict[(del_str, match_str)] = 0.05
        new_dict[(next_match, match_str)] = 0.9

        # transition from insert state I
        new_dict[(insert_str, insert_str)] = 0.5  # self-cycle
        new_dict[(next_match, insert_str)] = 0.5

        # transition from delete state
        if m_index < 139:
            new_dict[(next_del, del_str)] = 0.5
            new_dict[(two_ahead_match, del_str)] = 0.5
        else:  # no del->next_del transition
            # no del->match transition
            # only transition: d139 -> END
            new_dict[("END", del_str)] = 1.0

    # transitions from M140. Last match state
    # no delete, no (Match, Match)
    # transitions: (END, M140), (I140, M140)
    # transition from match state M
    new_dict[("Ir140", "Mr140")] = 0.05
    new_dict[("END", "Mr140")] = 0.95
    # I140 transitions
    new_dict[("Ir140", "Ir140")] = 0.5
    new_dict[("END", "Ir140")] = 0.5

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