#!/usr/bin/env python
#generates transition probs to and from delete, insert, match (sequence), and RNA_insert positions
#TODO MAKE PROB OF RNA SELF-CYCLE DEPENDENT ON AVG LENGTH
#for now, just a high prob: 0.85

#FORMAT:  dictionary of tuples: {(A,B):0.05, (C,E):0.5}
#in this example, we have prob of A given B is 0.05, prob of C given E is 0.5
import numpy as np

avg_lengths = [1641.4098732872096, 1539.26103329929, 1652.740349292529]
cells = ['2_B01', '3_C01', '4_D01']

#only difference in cell lines is prob of self-cycle-of  and from  RNA_insert
transition_dict_list = []#entry 0 is dict for 2B01, 1 for 3C01, 2 for 4D01
for i in np.arange(3):
    cell = cells[i]
    new_dict = {}
    #transition from match state M
    new_dict[("I","M")] = 0.05
    new_dict[("D","M")] = 0.05
    new_dict[("M","M")] = 0.9
    #transition from insert state I
    new_dict[("I","I")] = 0.5
    new_dict[("M","I")] = 0.5
    #transition from delete state D
    new_dict[("D","D")] = 0.5
    new_dict[("M","D")] = 0.5
    #transition from RNA-insert state R
    p = 0.85
    new_dict[("R","R")] = p
    new_dict[("","R")] = 1 - p - 0.05
    new_dict[("D","R")] = 0.05
    transition_dict_list.append(new_dict)


