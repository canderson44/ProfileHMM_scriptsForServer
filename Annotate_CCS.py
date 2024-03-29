#!/usr/bin/env python
# Purpose: pipeline so that given a list of CCSs, for each label find where the 5p barcode, 3p barcode, and adapter lie 

import CCS_Filter as cf
import numpy as np
import random as rd

SCORE_SIG_THRESHOLD = 0.9999999999 #score sig must be > this to be kept


#adaptation of CCS_Filter's score alignment, except this time: 
#RETURNS list of score, tuple of start and end coordinates, and string traceback alignment (just reference portion)
#FOR NOW: linear gap penalty d
#parameters: gap penalty d, barcode/reference(vertical), seq y(horizontal),
    #score matrix, reference character
#reference character: when returning alignment, replace all characters in seq y
    # with the reference character. Then only return the line of the alignment for 
    # the reference
    #example of return: 55555--5----------------
    #if the reference is the 5p barcode
#initialization: 0s in first row; (i,0) = -i*d
#Traceback favors iy>m>ix becuase more likely to have gaps in ref than 
    #y since reference shorter than y
def glocal_alignment(gap_penalty, sequence, reference, score_dict, reference_char):
    gap_penalty = np.abs(gap_penalty)
    
    #matrix
    num_rows = len(reference) +1
    num_cols = len(sequence) +1
    F = np.zeros((num_rows, num_cols))
    T = np.zeros((num_rows,num_cols)) #traceback
    
    #labels for tracing, prioritizing iy>m>ix
    IY, M, IX = 3, 2, 1
    T_list = [IY, M, IX]
    
    #initialization",
    #row 0 already initialized to zeros
    #keep T[0][j] as zeros because can stop anywhere there
    for i in np.arange(1,num_rows):
        F[i][0] = -1*i*gap_penalty
        T[i][0] = IX
        
    for i in np.arange(1,num_rows): # i represents the row
        for j in np.arange(1,num_cols): #j represents the col
            m = F[i-1][j-1] + score_dict[(reference[i-1],sequence[j-1])]
            ix = F[i-1][j] - gap_penalty
            iy = F[i][j-1] - gap_penalty
            F[i][j], T[i][j] = max(list(zip([iy,m,ix],T_list)))
    #now matrix is filled"
    last_row_array = np.array(F[num_rows-1])

    #max_indices = np.where(last_row_array == max(F[num_rows-1])) #only use highest score
    max_indices_list = []
    score_list = []
    # for elem in max_indices:
    #     for entry in elem:
    #         max_indices_list.append(entry)

    #populate max_indices_list with indices where scoresig > threshold
    for i in np.arange(len(F[num_rows-1])):
        score = F[num_rows-1][i]
        scoreSig = cf.calc_score_sig(score)
        if scoreSig > SCORE_SIG_THRESHOLD:
            max_indices_list.append(i)
    #sort max_indices_list with highest scoring indices first
    max_indices_list.sort(key=lambda index:F[num_rows-1][index], reverse=True)
    for index in max_indices_list:
        #print("index is", index)
        #print("adding to score_list this score: ", F[num_rows-1][index])
        score_list.append(F[num_rows-1][index])
    # print("number of columns:", len(F[num_rows-1]))
    # print("last row of F", F[num_rows-1])
    # print("max_indices_list", max_indices_list)
    # print("score_list", score_list)
    x_alignments = []
    y_alignments = []
    coords_list = [] #list of start indices in sequence string
    for trace_start_col in max_indices_list:
        #print("trace_start_col is ",trace_start_col)
        #start traceback
        x_align = [] #reference
        y_align = [] #seq
        seq_length = len(sequence)
        y = sequence
        x = [reference_char]*len(reference)
        i_fin, j_fin = num_rows-1, trace_start_col
        i = i_fin
        j = j_fin
        traceback = T[i][j]
        current = ''
        if j_fin < seq_length: #everything after j_fin is gapped for
            for n in np.arange(seq_length - j_fin):
                x_align.append('-')
                y_align.append(y[seq_length-n-1])
        if traceback == M:
            i -= 1
            j -= 1
            y_align.append(y[j])
            x_align.append(x[i])
        elif current == IY:
            j -= 1
            x_align.append('-')
            y_align.append(y[j])
        else: #current == IX
            i -= 1
            x_align.append(x[i])
            y_align.append('-')

        current = T[i][j]
        initial_coord = -1
        while i>0 and j>0:
            if current == M:
                i -= 1
                j -= 1
                y_align.append(y[j])
                x_align.append(x[i])
            elif current == IY:
                j -= 1
                x_align.append('-')
                y_align.append(y[j])
            else: #current == IX
                i -= 1
                x_align.append(x[i])
                y_align.append('-')
            current = T[i][j]

        if j >0:
            while j >0: #gapping in beg of sequence; region starts before sequence
                j -= 1
                x_align.append('-')
                y_align.append(y[j])
        x_align_copy = x_align.copy() #copy of alignment in reverse order
        x_align.reverse()
        y_align.reverse()
        #help with troubleshooting coords
        # x_str = "".join(x_align)
        # alt_start = x_str.find(reference_char) #first index of occurence
        # alt_end_0 = len(x_align) - 1 - x_align[::-1].index(reference_char)
        # for space in np.arange(len(x_align)):
        #     if x_align[space] == '-':
        #         x_align[space] = str(space)
        #print("region alignment is", x_align)
        start_coord = x_align.index(reference_char)
        stop_coord = len(x_align) - 1 - x_align_copy.index(reference_char)
        if (start_coord, stop_coord, reference_char) not in coords_list:
            #print("appending coords", (start_coord, stop_coord, reference_char))
            # print()
            coords_list.append((start_coord,stop_coord, reference_char))
            x_alignments.append(x_align)
            y_alignments.append(y_align)
        else:
            pass
            #print("coords " + str((start_coord, stop_coord, reference_char)) + "already in coords list. don't add this alignment.")

        
    return_list = list(zip(score_list,coords_list,x_alignments))
    return return_list

#important strings
nucleotides = ['A','C','G','T']
barcodes_list = []
with open('/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/barcodes_for_profileHMM.fasta') as f:
    for line in f:
        line = line.strip()
        if line[0]!='>':
            barcodes_list.append(line)
barcode4 = barcodes_list[1]
barcode3 = barcodes_list[2]
barcode2 = barcodes_list[3]
fivePBarcode = barcodes_list[0]
fivePBarcodeRC = cf.gen_rev_complement(fivePBarcode)
barcode2RC = cf.gen_rev_complement(barcode2)
barcode3RC = cf.gen_rev_complement(barcode3)
barcode4RC = cf.gen_rev_complement(barcode4)
#adapter
two_adapter = ''
three_adapter = ''
four_adapter = ''
filename2 = '/tier2/deweylab/data/thomson_lab/Pacbio/2_B01/m54178_180915_120213.adapters.fasta'
with open(filename2) as f: 
    for line in f: 
        line = line.strip()
        if line[0]!='>':
            two_adapter = line

filename3 = '/tier2/deweylab/data/thomson_lab/Pacbio/3_C01/m54178_180916_082219.adapters.fasta'
with open(filename3) as f: 
    for line in f: 
        line = line.strip()
        if line[0]!='>':
            three_adapter = line

filename4 = '/tier2/deweylab/data/thomson_lab/Pacbio/4_D01/m54178_180917_044252.adapters.fasta'            
with open(filename4) as f: 
    for line in f: 
        line = line.strip()
        if line[0]!='>':
            four_adapter = line
two_adapter_RC = cf.gen_rev_complement(two_adapter)
three_adapter_RC = cf.gen_rev_complement(three_adapter)
four_adapter_RC = cf.gen_rev_complement(four_adapter)

#mutations
twoAdapter_mutations = cf.gen_mutated_barcodes(two_adapter)
threeAdapter_mutations = cf.gen_mutated_barcodes(three_adapter)
fourAdapter_mutations = cf.gen_mutated_barcodes(four_adapter)
bar2_mutations = cf.gen_mutated_barcodes(barcode2)
bar2RC_mutations = cf.gen_mutated_barcodes(barcode2RC)
bar3_mutations = cf.gen_mutated_barcodes(barcode3)
bar3RC_mutations = cf.gen_mutated_barcodes(barcode3RC)
bar4_mutations = cf.gen_mutated_barcodes(barcode4)
bar4RC_mutations = cf.gen_mutated_barcodes(barcode4RC)
fivePbar_mutations = cf.gen_mutated_barcodes(fivePBarcode)
fivePbarRC_mutations = cf.gen_mutated_barcodes(fivePBarcodeRC)

#generates list of test strings given barcode, adapter, and respective rev compliments
#PARAMETERS: barcode_index: index of 3' barcode:
                # if cell is 2_B01, index is 2; 3_C01 index is 3; 4_D01 index is 4
def gen_test_strings(barcode_index, barcode, barcodeRC, adapter, adapterRC):
    if barcode_index == 2:
        # list of substituted, inserted, deleted mutations on barcode
        threeBar_mutations = bar2_mutations
        threeBarRC_mutations = bar2RC_mutations
        adapter_mutations = twoAdapter_mutations
    elif barcode_index == 3:
        # list of substituted, inserted, deleted mutations on barcode
        threeBar_mutations = bar3_mutations
        threeBarRC_mutations = bar3RC_mutations
        adapter_mutations = threeAdapter_mutations
    else: #barcode_index == 4
        # list of substituted, inserted, deleted mutations on barcode
        threeBar_mutations = bar4_mutations
        threeBarRC_mutations = bar4RC_mutations
        adapter_mutations = fourAdapter_mutations
    test_strings = []
    #Ar3r5r
    test_ar3r5r = adapterRC + barcodeRC + fivePBarcodeRC
    test_strings.append(test_ar3r5r)

    #A53A
    test_a53a = adapter + "".join(rd.choices(nucleotides, k=50))+ fivePBarcode+\
                "".join(rd.choices(nucleotides, k=100)) + barcode + \
                "".join(rd.choices(nucleotides, k=10)) + adapter
    test_strings.append(test_a53a)
    #5rA3Ar
    test_5ra3ar = fivePBarcodeRC + adapter + "".join(rd.choices(nucleotides, k=10)) +\
                  barcode +"".join(rd.choices(nucleotides, k=100)) +adapterRC
    test_strings.append(test_5ra3ar)
    #start5
    test_start5= fivePBarcode + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_start5)
    #start5Sub
    test_start5Sub =  fivePbar_mutations[0]+ "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_start5Sub)
    #start5Insert
    test_start5Insert = fivePbar_mutations[1] + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_start5Insert)
    #start5Delete
    test_start5Del = fivePbar_mutations[2] + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_start5Del)
    #A_rand_A
    test_a_rand_a = adapter + "".join(rd.choices(nucleotides, k=900)) + adapter
    test_strings.append(test_a_rand_a)
    #Astartsub
    test_astartSub = adapter_mutations[0] + "".join(rd.choices(nucleotides, k=900))
    test_strings.append(test_astartSub)
    #Ainsert
    test_astartInsert = adapter_mutations[1] + "".join(rd.choices(nucleotides, k=900))
    test_strings.append(test_astartInsert)
    #Adelete
    test_astartDel = adapter_mutations[2] + "".join(rd.choices(nucleotides, k=900))
    test_strings.append(test_astartDel)
    #rand
    test_random="".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_random)
    #end in 3
    test_endIn3pBarcode= "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_endIn3pBarcode)

    # following check that keep all above threshold, not just highest scoring for some region
    # also, since beginning of 3 is 5' barcode exactly, makes sure 5' eliminated, not intended 3'

    #just 3 sub
    test_strings.append(threeBar_mutations[0])

    #3 sub 3 sub
    test_strings.append(threeBar_mutations[0] + "".join(rd.choices(nucleotides, k=1000)) + threeBar_mutations[0])
    #3sub 3
    test_strings.append(threeBar_mutations[0] + "".join(rd.choices(nucleotides, k=10)) + barcode)
    #3insert 3
    test_strings.append(threeBar_mutations[1]+ "".join(rd.choices(nucleotides, k=100)) + barcode)
    #3del 3
    test_strings.append(threeBar_mutations[2]+ "".join(rd.choices(nucleotides, k=100)) + barcode)

    #A_3_3
    test_a33 = test_aSub3_3 = adapter + "".join(rd.choices(nucleotides, k=100)) + barcode + "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_a33)

    #A_Sub3_3
    test_aSub3_3 = test_aSub3_3 = adapter + "".join(rd.choices(nucleotides, k=100)) + threeBar_mutations[0] + "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_aSub3_3)
    #A_Insert3_3
    test_aInsert3_3 = adapter + "".join(rd.choices(nucleotides, k=100)) + threeBar_mutations[1] + "".join(
        rd.choices(nucleotides, k=100)) + barcode
    test_strings.append(test_aInsert3_3)
    #A_Del3_3
    test_aDel3_3 = adapter + "".join(rd.choices(nucleotides, k=100)) + threeBar_mutations[2] + "".join(
        rd.choices(nucleotides, k=100)) + barcode
    test_strings.append(test_aDel3_3)
    return test_strings

#given a sequence, and a list of lists of [reference, reference character], 
#annotates the sequence with the location(s) of the reference(s)
#does references and their reverse complements, with the RCs having "r" 
    #appended to the reference character
#only includes reference(s) with a score significance of >0.99999999

#RETURN: of references that pass score significance filter, returns:

      #case 1: justCoords == True
            #list of tuples: each tuple is (start_coord, stop_coord, reference_char)
            #where start and stop coords are start and stop indices within sequence string where region of interest aligns

       #case 2: justCoords == False
 
            #returns fusion of the alignment strings returned by glocal alignment function 

#NOTE: of a region and its reverse complement, will include a max of one of them, whichever is higher (but only if they pass the score significance fitler) 

#PARAMETERS
    #sequence: string
    #ref_list: list of lists
        #each entry is [reference(string), reference character]
def annotate_seq(sequence, ref_list, justCoords=False):

    final_coord_list = [] #for case 1
    fused_list = ['-'] * len(sequence) #for case 2
    sorted_annotations_list = [] #list of annotated regions sorted by scoreSig
        #used for our greedy selection of regions to keep
    selected_annotations_list = []  # holds annotations we greedily choose to keep

    for refPair in ref_list:
        #print("ref_pair is" ,refPair)
        ref = refPair[0] #region of interest
        ref_char = refPair[1]
        #print("this ref is: ", ref_char)
        refRC = cf.gen_rev_complement(ref)
        refRC_char = ref_char + 'r'

        #print("in annotate_seq. ref_char is", ref_char)
        ########
        ########
        ########
        #first: forward
        ########
        ########
        ########

        #gloc is tuple of lists: (score_list,coords_list,x_alignments)
        gloc = glocal_alignment(gap_penalty=1, sequence=sequence,
                                    reference=ref, score_dict = cf.score_dict,
                                     reference_char = ref_char)
        #print("glocal alignment tuple is", gloc)
        #this_ref_list: list of tuples: [(score,alignment string), ...]
        this_ref_list = [(gloc[n][0],gloc[n][2]) for n in np.arange(len(gloc))]

        #current_coords: list of tuples: [(start,stop,char),...]
        current_coords = [gloc[n][1] for n in np.arange(len(gloc))]
        ref_scores = []
        ref_alignments = []
        #for now, keep all
        for pair in this_ref_list:
            ref_scores.append(pair[0])
            ref_alignments.append(pair[1])

        ########
        ########
        ########
        # now: reverse complement
        ########
        ########
        ########
        rc_gloc = glocal_alignment(gap_penalty=1, sequence=sequence,
                                              reference=refRC, score_dict = cf.score_dict, reference_char = refRC_char)
        rc_list = [(rc_gloc[n][0],rc_gloc[n][2]) for n in np.arange(len(rc_gloc))] #score, alignment
        rc_current_coords = [rc_gloc[n][1] for n in np.arange(len(rc_gloc))]
        rc_scores = []
        rc_alignments = []

        #for now, keep all in rc_list
        for pair in rc_list:
            rc_scores.append(pair[0])
            rc_alignments.append(pair[1])

        ########
        ########
        ########
        # calc score sigs, determine which to use: ref or rc
        ########
        ########
        ########
        zeroRef = False
        zeroRC = False
        useRef = False
        allEmpty = False
        if len(ref_scores) > 0:
        #    print("ref scores list is", ref_scores)
            ref_scoreSig_list = [cf.calc_score_sig(ref_score) for ref_score in ref_scores]
        #    print("ref scoreSig list is ", ref_scoreSig_list)
            max_ref_scoresig = max(ref_scoreSig_list)
        else:
            zeroRef = True

        if len(rc_scores) >0:
            rc_scoreSig_list = [cf.calc_score_sig(rc_score) for rc_score in rc_scores]
            max_rc_scoresig = max(rc_scoreSig_list)

        else:
            zeroRC = True
        if zeroRC == False and zeroRef == False:
            #print("zeroRC and zeroRef are both false")
            if max_ref_scoresig > max_rc_scoresig:
                useRef = True
        elif zeroRC == False and zeroRef == True:
            #print("zeroRC si false, zeroRef is true")
            useRef = False
        elif zeroRC == True and zeroRef == False:
            #print("zeroRC is true, zeroRef is false")
            useRef = True
        else: #both are true
            #print("allEmpty is true")
            allEmpty = True


        #now: case of all empty (no useable coords) vs not empty
        if allEmpty == False:
            #print("in allEmpty=False case")
            ########
            ########
            ########
            # Add annotated regions to list
            ########
            ########
            ########
            # conditions: 1. only add to list if score_sig above threshold and its
            #                   appropriate direction (forward or rc)
            #             2. format: (scoreSig, coord tuple, alignment string)
            #             3. sort based on scoreSig
            #print("for this region type, before filtering by score sig, list of potential regions is", current_coords)
            if useRef == True:
                for index in np.arange(len(ref_scoreSig_list)):
                    #print("this potential annotation is: ", current_coords[index])
                    #print("this score sig is: ", ref_scoreSig_list[index])
                    if ref_scoreSig_list[index]  > SCORE_SIG_THRESHOLD:
                        #print("adding this region")
                        sorted_annotations_list.append((ref_scoreSig_list[index],
                                        current_coords[index], ref_alignments[index]))
            else: #use rc
                for rc_index in np.arange(len(rc_scoreSig_list)):
                    #print("this potential annotation is: ", rc_current_coords[rc_index])
                    #print("this score sig is: ", ref_scoreSig_list[rc_index])
                    if rc_scoreSig_list[rc_index] > SCORE_SIG_THRESHOLD:
                        sorted_annotations_list.append((rc_scoreSig_list[rc_index],
                                rc_current_coords[rc_index], rc_alignments[rc_index]))
    #print("before sorting, sorted annotations list is")
    # for elem in sorted_annotations_list:
    #     print_str = "(" + str(elem[0]) + str(elem[1]) + ")"
    #     print(print_str)
    if len(sorted_annotations_list) > 0:
        #print("pre-sorted sorted list is not empty")
        ########
        ########
        ########
        # Sort the list
        ########
        ########
        ########
        # now we'll actually sort the sorted_annotations_list
        # conditions: sort based on 0th element of each list element. Sort in descending order
        #           i.e. sorted by scoreSig, highest scoreSig is first in list
        sorted_annotations_list.sort(key=lambda elem: elem[0], reverse=True)
        # print("FINAL SORTED-ANNOTATIONS-LIST")
        # for elem in sorted_annotations_list:
        #     print_str = "(" + str(elem[0]) + str(elem[1]) + ")"
        #     print(print_str)
        # print()
        ########
        ########
        ########
        # Greedily select from now-sorted sorted_annotations_list
        ########
        ########
        ########
        #######
        # conditions for selection:
        #     1. select and remove from top of list (highest scoreSig)
        #     2. add to list of selected if new annotation's coords don't overlap with
        #               any of the coords already in teh list of selected


        # recall: elements of sorted_annotations_list are of form
        #                           (scoreSig, coord tuple, alignment string)
        # selected_annotations_list has same format

        #first: pop first element into the selected list, so we have something to compare to
        if len(sorted_annotations_list) >0:
            selected_annotations_list.append(sorted_annotations_list.pop(0))
            #print("first added to selected annotations: ", selected_annotations_list[0][1])

        # stop once we've considered all viable annotations
        while(len(sorted_annotations_list) > 0):
            maybe_annotation_tuple = sorted_annotations_list.pop(0)
            #print("maybe annotation: ", maybe_annotation_tuple[1])
            maybe_coords = maybe_annotation_tuple[1]
            maybe_start = maybe_coords[0]
            maybe_stop = maybe_coords[1]
            addMaybe = True #if still true by end of for loop below, we'll add the maybe tuple
            for selected_tuple in selected_annotations_list:
                selected_coords = selected_tuple[1]
                selected_start = int(selected_coords[0])
                selected_stop = int(selected_coords[1])
                #let's check for overlap
                #note: +1 in the np.arange stop because stop coords are inclusive
                maybe_range = [n for n in np.arange(maybe_start, maybe_stop + 1)]
                selected_range = [n for n in np.arange(selected_start, selected_stop + 1)]
                #only want to keep maybe if no overlap between selected and maybe coords
                if ((selected_start in maybe_range) or (selected_stop in maybe_range) or
                        (maybe_start in selected_range) or (maybe_stop in selected_range)):
                    #print("this maybe won't be added")
                    addMaybe = False
            if(addMaybe):
                #print("this maybe is added")
                # then we do want to keep the maybe_annotation
                selected_annotations_list.append(maybe_annotation_tuple)

    # print("SELECTED_ANNOTATIONS_LIST")
    # for elem in selected_annotations_list:
    #     print_str = "(" + str(elem[0]) + str(elem[1]) + ")"
    #     print(print_str)
    # print()
    ########
    ########
    ########
    # wrap up coords or annotation strings to return
    ########
    ########
    ########
    #######
    # here we split into case 1 and case 2

    #case 1
    if justCoords:
       #want to return list of tuples: [(start,stop,ref_char),....]
       # selected tuple format: (scoreSig, coord tuple, alignment string)
       # cord tuple format: (start,stop,char)
        for selected_tuple in selected_annotations_list:
           this_coord_pair = selected_tuple[1]
           final_coord_list.append(this_coord_pair)
        #print("returning final coords list:", final_coord_list)
        return final_coord_list

    #case 2
    else: #justCoords == False
        #want to return a string: fuse(fused_list,alignments[i])
        # selected tuple format: (scoreSig, coord tuple, alignment string)
        for selected_tuple in selected_annotations_list:
            fuse(fused_list, selected_tuple[2])
        return "".join(fused_list)

    
#given current fused list and list of what to fuse, fuses the "to_fuse" contents to the fused list
def fuse(fused, to_fuse):
    len_fused = len(fused)
    for i in np.arange(len_fused):
        if to_fuse[i] != '-':
            if fused[i] != '-':
                fused[i] = '!'
            else:
                fused[i] = to_fuse[i]

#test it!
#test fuse:
def test_fuse(filler):
    filler = '-'
    fused = [filler] * 5
    toFuse0 = [filler,filler,filler,'0','0']
    fuse(fused,toFuse0)
    toFuse1 = [filler,'1','1','1',filler]
    fuse(fused, toFuse1)
    assert fused == [filler, '1', '1', '1', '0']

#test annotate_seq
ref_list_2 = [[fivePBarcode,'5'],[barcode2,'2'], [two_adapter,'A'] ]
ref_list_3 = [[barcode3,'3'], [three_adapter,'A'], [fivePBarcode,'5']]
ref_list_4 = [[barcode4,'4'], [four_adapter,'A'], [fivePBarcode,'5']]
def test_annotateseq(ref_list, seq_list):
    for seq in seq_list:
        print("CALC FINAL COORDS")
        final_coords = annotate_seq(seq, ref_list, justCoords=True)
        print("CALC FINAL FUSED STRING")
        final_fused_list = annotate_seq(seq,ref_list,justCoords=False)
        print("format: final_coord_list: [(start,stop,char),(start,stop,char),etc]")
        print(final_coords)
        print("format: fused string")
        print(final_fused_list)
        print("################\n##############\n##############\n#############\n\n")





# test annotateseq:
#test_strings_2 = gen_test_strings(2,barcode2, barcode2RC, two_adapter, two_adapter_RC)
#test_strings_3 = gen_test_strings(3,barcode3, barcode3RC, three_adapter, three_adapter_RC)
#test_strings_4 = gen_test_strings(4,barcode4, barcode4RC, four_adapter, four_adapter_RC)
#test_annotateseq(ref_list_2,test_strings_2)
#test_annotateseq(ref_list_3,test_strings_3)
#test_annotateseq(ref_list_4, test_strings_4)
