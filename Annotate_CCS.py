#!/usr/bin/env python
# Purpose: pipeline so that given a list of CCSs, for each label find where the 5p barcode, 3p barcode, and adapter lie 

import CCS_Filter as cf
import numpy as np
import random as rd

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
    max_indices = np.where(last_row_array == max(F[num_rows-1]))
    max_indices_list = []
    for elem in max_indices:
        for entry in elem:
            max_indices_list.append(entry)
    score_list = [F[num_rows-1][i] for i in max_indices_list]
    x_alignments = []
    y_alignments = []
    coords_list = [] #list of start indices in sequence string
    for trace_start_col in max_indices_list:
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
        if j_fin < seq_length:
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
            while j >0:
                j -= 1
                x_align.append('-')
                y_align.append(y[j])

        x_align.reverse()
        y_align.reverse()
        x_str = "".join(x_align)
        start_coord = x_str.find(reference_char) #first index of occurence
        stop_coord = x_str.rfind(reference_char)#last index of occurence (searches backwards)
        coords_list.append((start_coord,stop_coord))
        x_alignments.append(x_align)
        y_alignments.append(y_align)
        
    return_list = list(zip(score_list,coords_list,x_alignments))
    return return_list

#important strings
nucleotides = cf.nucleotides
fivePBarcode = cf.barcodes_list[0]
fivePBarcodeRC = cf.gen_rev_complement(fivePBarcode)
barcode2 = cf.barcode_2
barcode3 = cf.barcode_3
barcode4 = cf.barcode_4
barcode2RC = cf.barcode_2_rc
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

#generates list of test strings given barcode, adapter, and respective rev compliments
def gen_test_strings(barcode, barcodeRC, adapter, adapterRC):    
    test_strings = []
    #First: exact matches
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
                  barcode+"".join(rd.choices(nucleotides, k=100)) +adapterRC
    test_strings.append(test_5ra3ar)
    #start5
    test_start5= fivePBarcode + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_start5)
    #A_rand_A
    test_a_rand_a = adapter + "".join(rd.choices(nucleotides, k=900)) + adapter
    test_strings.append(test_a_rand_a)
    #rand
    test_random="".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_random)
    #end in 3
    test_endIn3pBarcode= "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_endIn3pBarcode)
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
    for refPair in ref_list:
#        print("ref_pair is" ,refPair)
        ref = refPair[0]
        ref_char = refPair[1]
        refRC = cf.gen_rev_complement(ref)
        refRC_char = ref_char + 'r'
        gloc = glocal_alignment(gap_penalty=1, sequence=sequence,
                                              reference=ref, score_dict = cf.score_dict,
                                     reference_char = ref_char)
#        print("len of gloc is: ", len(gloc))
#        print("elements of gloc: ", gloc)
        this_ref_list = [(gloc[n][0],gloc[n][2]) for n in np.arange(len(gloc))]
        current_coords = gloc[1]
        ref_scores = []
        ref_alignments = []
        for pair in this_ref_list:
            ref_scores.append(pair[0])
            ref_alignments.append(pair[1])
        #now reverse complement
        rc_gloc = glocal_alignment(gap_penalty=1, sequence=sequence,
                                              reference=refRC, score_dict = cf.score_dict, reference_char = refRC_char)
        rc_list = [(rc_gloc[n][0],rc_gloc[n][2]) for n in np.arange(len(rc_gloc))] #score, alignment
        rc_current_coords = [rc_gloc[n][1] for n in np.arange(len(rc_gloc))]
        rc_scores = []
        rc_alignments = []
        for pair in rc_list:
            rc_scores.append(pair[0])
            rc_alignments.append(pair[1])

        ref_scoreSig_list = [cf.calc_score_sig(ref_score) for ref_score in ref_scores]
        rc_scoreSig_list = [cf.calc_score_sig(rc_score) for rc_score in rc_scores]
        useRef = False
        for i in np.arange(len(ref_scoreSig_list)):
            for j in np.arange(len(rc_scoreSig_list)):
                if ref_scoreSig_list[i] > rc_scoreSig_list[j]:
                    useRef = True
                else:
                    useRef = False
        # here we split into case 1 and case 2
        if justCoords:
            if useRef == True:
                for i in np.arange(len(ref_scoreSig_list)):
                    if ref_scoreSig_list[i]>0.999999999:
                        final_coord_list.append((current_coords[0],current_coords[1],ref_char))
            else: #useRef == False
                for i in np.arange(len(rc_scoreSig_list)):
                    if rc_scoreSig_list[i]>0.999999999:
                        final_coord_list.append((rc_current_coords[0], rc_current_coords[1], refRC_char))

            return final_coord_list

        else: #justCoords == False; casae 2
            if useRef == True:
                for i in np.arange(len(ref_scoreSig_list)):
                    if ref_scoreSig_list[i]>0.999999999:
                        fuse(fused_list,ref_alignments[i])
            else: #useRef == False
                for i in np.arange(len(rc_scoreSig_list)):
                    if rc_scoreSig_list[i]>0.999999999:
                        fuse(fused_list,rc_alignments[i])
            return_str = "".join(fused_list)
            return return_str
    
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
        print(annotate_seq(seq, ref_list, justCoords=True))
        print(annotate_seq(seq,ref_list,justCoords=False))
        print()





# test annotateseq:
test_strings_2 = gen_test_strings(barcode2, barcode2RC, two_adapter, two_adapter_RC)
test_strings_3 = gen_test_strings(barcode3, barcode3RC, three_adapter, three_adapter_RC)
test_strings_4 = gen_test_strings(barcode4, barcode4RC, four_adapter, four_adapter_RC)
test_annotateseq(ref_list_2,test_strings_2)
# test_annotateseq(ref_list_3,test_strings_3)
# test_annotateseq(ref_list_4, test_strings_4)
