#!/usr/bin/env python

# Program for separating useable ZMW CCSs from unusable. Useable if don't begin with 3' barcode and do end with 3' barcode. 
# 
# kept sequence must have 1 and 2 as follows:
#     1. 
#         3' barcode in end region of sample (sample is forward transcript)
#         OR
#         3' barcode reverse complement in beginning region of sample (sample is rev compliment)
#     AND
#     2. 
#         no 3' barcode in beginning region of sample (forward transcript)
#         no 3' barcode reverse complement in end region of sample (rc transcript)
import traceback
import os
import numpy as np
import random as rd
from time import sleep

#given DNA nucleotide, returns its complement
def give_nuc_complement (nucleotide):
    if nucleotide == "A":
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T': 
        return 'A'
    else: #
        return '?'


#generates reverse compliment. Returns rc in standard 5' to 3' direction
def gen_rev_complement(sequence):
    complement_list = []
    for nuc in sequence:
        complement_list.append(give_nuc_complement(nuc))
    #have reverse complement, but it's backwards
    complement_list.reverse()
    return "".join(complement_list)


#gets your beginning and ending sample from sequence. Either one is 1.5 * passed length. 
def get_beg_and_end(sample, barcode_length):
    segment_length = int(1.5*barcode_length)
    beg = sample[:segment_length]
    end = sample[len(sample)-segment_length:]
    return [beg,end]

barcodes_list = []
with open('/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/barcodes_for_profileHMM.fasta') as f: 
    for line in f:
        line = line.strip()
        if line[0]!='>':
            barcodes_list.append(line)
barcode_4 = barcodes_list[1]
barcode_3 = barcodes_list[2]
barcode_2 = barcodes_list[3]
barcode_2_rc = gen_rev_complement(barcode_2)

#adapter
two_adapter = ''
with open('/tier2/deweylab/data/thomson_lab/Pacbio/2_B01/m54178_180915_120213.adapters.fasta') as f:
    for line in f: 
        line = line.strip()
        if line[0]!='>':
            two_adapter = line
            

#test "CCSs"
nucleotides = ['A','C','G','T']

def gen_mutated_barcodes(barcode):
    barcode_list = list(barcode)
    rand_subs_list = rd.choices(np.arange(len(barcode)), k=int(len(barcode) / 10))  # pacbio MUCH more accurate
    # here, error in sequencing is 10%
    # pacbio has error of about 0.001%
    # source: https://www.pacb.com/uncategorized/a-closer-look-at-accuracy-in-pacbio/
    sub_list = barcode_list.copy()
    rand_nuc_list = rd.choices(nucleotides, k=len(rand_subs_list))

    for i in np.arange(len(rand_nuc_list)):
        sub_list[i] = rand_nuc_list[i]
    sub_str = "".join(sub_list)

    insert_list = barcode_list.copy()
    for i in np.arange(len(rand_nuc_list)):
        insert_list.insert(i, rand_nuc_list[i])
    insert_str = "".join(insert_list)

    delet_list = barcode_list.copy()
    for i in np.arange(int(len(rand_nuc_list))):
        del barcode_list[i]
    delet_str = "".join(delet_list)
    return [sub_str, insert_str, delet_str]

bar2_mutations = gen_mutated_barcodes(barcode_2)
bar2RC_mutations = gen_mutated_barcodes(barcode_2_rc)
bar3_mutations = gen_mutated_barcodes(barcode_3)
bar3RC_mutations = gen_mutated_barcodes(gen_rev_complement(barcode_3))
bar4_mutations = gen_mutated_barcodes(barcode_4)
bar4RC_mutations = gen_mutated_barcodes(gen_rev_complement(barcode_4))
fivePbar_mutations = gen_mutated_barcodes(barcodes_list[0])
fivePbarRC_mutations = gen_mutated_barcodes(gen_rev_complement(barcodes_list[0]))

#generates a list of test strings given 3p barcode and its rev complement
#list is: 
#exact matches to barcode:
##endsIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3
#substituted, inserted, and deleted versions of barcode:
#endIn3RC * 3, begIn3RC * 3, endIn3 * 3, begIn3 * 3, endIn3 * 3
def gen_test_strings(barcode, barcode_rc, threePbar_mutations, threePbarRC_mutations):
    test_strings = []

    # #First: exact matches
    # #exclude
    # test_end_in_3pBarcodeRC = "".join(rd.choices(nucleotides, k=1000)) + barcode_rc
    # test_strings.append(test_end_in_3pBarcodeRC)
    # #keep
    # test_beg_3pBarcodeRC = barcode_rc + "".join(rd.choices(nucleotides, k=1000))
    # test_strings.append(test_beg_3pBarcodeRC)
    # #keep
    # test_5pBarcodeRandom3pBarcode = barcodes_list[0] + "".join(rd.choices(nucleotides, k=1000)) + barcode
    # test_strings.append(test_5pBarcodeRandom3pBarcode)
    # #exclude
    # test_3pBarcodeThenRandom= barcode + "".join(rd.choices(nucleotides, k=(1000)))
    # test_strings.append(test_3pBarcodeThenRandom)
    # #exclude
    # test_random="".join(rd.choices(nucleotides, k=(1000)))
    # test_strings.append(test_random)
    # #keep
    # test_endIn3pBarcode= "".join(rd.choices(nucleotides, k=1000)) + barcode
    # test_strings.append(test_endIn3pBarcode)


    #Next: mutated matches



    #exclude
    test_end3RC_sub = "".join(rd.choices(nucleotides, k=1000)) + threePbarRC_mutations[0]
    test_end3RC_insert = "".join(rd.choices(nucleotides, k=1000)) + threePbarRC_mutations[1]
    test_end3RC_del = "".join(rd.choices(nucleotides, k=1000)) + threePbarRC_mutations[2]
    test_strings.append(test_end3RC_sub)
    test_strings.append(test_end3RC_insert)
    test_strings.append(test_end3RC_del)
    #keep
    test_beg3RC_sub = threePbarRC_mutations[0] + "".join(rd.choices(nucleotides, k=1000))
    test_beg3RC_insert = threePbarRC_mutations[1] + "".join(rd.choices(nucleotides, k=1000))
    test_beg3RC_del = threePbarRC_mutations[2] + "".join(rd.choices(nucleotides, k=1000))
    test_strings.append(test_beg3RC_sub)
    test_strings.append(test_beg3RC_insert)
    test_strings.append(test_beg3RC_del)

    #keep
    test_5pRand3p_sub = fivePbar_mutations[0] + "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[0]
    test_5pRand3p_insert = fivePbar_mutations[1] + "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[1]
    test_5pRand3p_del = fivePbar_mutations[2] + "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[2]
    test_strings.append(test_5pRand3p_sub)
    test_strings.append(test_5pRand3p_insert)
    test_strings.append(test_5pRand3p_del)

    #exclude
    test_beg3p_sub = threePbar_mutations[0] + "".join(rd.choices(nucleotides, k=1000))
    test_beg3p_insert = threePbar_mutations[1] + "".join(rd.choices(nucleotides, k=(1000)))
    test_beg3p_del = threePbar_mutations[2] + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_beg3p_sub)
    test_strings.append(test_beg3p_insert)
    test_strings.append(test_beg3p_del)
    
    #keep
    test_end3p_sub = "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[0]
    test_end3p_insert = "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[1]
    test_end3p_del = "".join(rd.choices(nucleotides, k=1000)) + threePbar_mutations[2]
    test_strings.append(test_end3p_sub) 
    test_strings.append(test_end3p_insert) 
    test_strings.append(test_end3p_del)
    
    return test_strings

test_2B01_strings = gen_test_strings(barcode_2, barcode_2_rc, bar2_mutations, bar2RC_mutations)
test_3C01_strings = gen_test_strings(barcode_3, gen_rev_complement(barcode_3), bar3_mutations, bar3RC_mutations)
test_4D01_strings = gen_test_strings(barcode_4, gen_rev_complement(barcode_4), bar4_mutations, bar4RC_mutations)

expected_keep_results = [False, True, True, False, False, True, False, False, False, True, True, True, 
                             True, True, True, False, False, False, 
                             True, True, True]

score_dict = {}
mismatch_penalty = -2
score_dict[('A','A')] = 1
score_dict[('A','C')] = mismatch_penalty
score_dict[('A','G')] = mismatch_penalty
score_dict[('A','T')] = mismatch_penalty
score_dict[('C','A')] = mismatch_penalty
score_dict[('C','C')] = 1
score_dict[('C','G')] = mismatch_penalty
score_dict[('C','T')] = mismatch_penalty
score_dict[('G','A')] = mismatch_penalty
score_dict[('G','C')] = mismatch_penalty
score_dict[('G','G')] = 1
score_dict[('G','T')] = mismatch_penalty
score_dict[('T','A')] = mismatch_penalty
score_dict[('T','C')] = mismatch_penalty
score_dict[('T','G')] = mismatch_penalty
score_dict[('T','T')] = 1


#performs SmithWaterman adaptation (glocal alignment). Returns score; no traceback
    #FOR NOW: linear gap penalty d
    #parameters: gap penalty, seq y (horizontal), barcode reference (vertical),  
    #initialization: 0s in first row; (i,0) = -i*d
    #sequence larger than barcode; barcode is reference. 
        #sequence is 1.5 * barcode length
def score_alignment(gap_penalty, sequence, barcode, score_dict):
    gap_penalty = np.abs(gap_penalty)
    #matrix
    num_rows = len(barcode) +1
    num_cols = len(sequence) +1
    F = np.zeros((num_rows, num_cols))
    #initialization",
    #row 0 already initialized to zeros
    for i in np.arange(1,num_rows):
        F[i][0] = -1*i*gap_penalty

    for i in np.arange(1,num_rows): # i represents the row
        for j in np.arange(1,num_cols): #j represents the col
            m = F[i-1][j-1] + score_dict[(barcode[i-1],sequence[j-1])]
            ix = F[i-1][j] - gap_penalty
            iy = F[i][j-1] - gap_penalty


            F[i][j] = max(m,ix,iy)
    #now matrix is filled\n",
    score = max(F[num_rows-1])
    return score

#compute s prime = score + log(pm/pr)
#score: alignment score
#pm: prior prob of m
#pr: prior prob of r
def compute_s_prime(score, pm, pr):
    return score + np.log(pm/pr)

#computes sigma function e^x / (1+e^x)
def compute_sigma(x):
    return np.exp(x)/(1.0+np.exp(x))


#calculates score significance with Bayesian approach: P(M|x,y)
#if P(M|x,y)>0.9, x and y related
#possibilities for sequence of interest: 3' adapter, 3' barcode, 
                        #RNA insert, 5' barcode, 5' adapter. Only one is M
#P(M) = 1/5 Match model
#P(R) = 4/5 Random model
def calc_score_sig(score):
    pr = 4.0/5
    pm = 1.0/5
    s_prime = compute_s_prime(score, pr, pm)  
    return compute_sigma(s_prime)



#given sequence, of CCS, 
#determines whether it should be kept
#True means keep, false means exclude

#kept sequence must have 1 and 2 as follows:
    #1. 
        #3' barcode in end region of sample (if sample is forward transcript)
        #OR
        #3' barcode reverse compliment in beginning region of sample (if sample is rev compliment)
    #AND
    #2. 
        #no 3' barcode in beginning region of sample (forward transcript)
        #no 3' barcode reverse complement in end region of sample (rc transcript)

#simplified. 
def should_we_keep(sample, barcode):
    segments = get_beg_and_end(sample, len(barcode))
    barcode_rc = gen_rev_complement(barcode)
    #beg
    #score_alignment(gap_penalty, sequence, barcode, score_dict):
    barcode_beg_score_sig = calc_score_sig(score_alignment(1,segments[0],barcode,score_dict))
    barcodeRC_beg_score_sig = calc_score_sig(score_alignment(1,segments[0],barcode_rc,score_dict))
    barcode_end_score_sig = calc_score_sig(score_alignment(1,segments[1],barcode,score_dict))
    
    if barcode_beg_score_sig >0.99999999:
        return False            
    elif  barcodeRC_beg_score_sig>0.99999999:
        return True
    #end                        
    elif barcode_end_score_sig > 0.99999999:
        return True
    else: #if we've gottent to this point, either the barcode isn't in teh 
        #sample and/or the barcode rev comp is at the end, which we don't want
        return False
        

#test it
#tests should_we_keep function given barcode and test strings
def test_should_we_keep(test_strings, barcode): 
    keep_list = []
    for entry in test_strings:
        keep_list.append(should_we_keep(entry, barcode))
#    print('results ',keep_list)
#    print('expected',expected_keep_results)
    
#test_should_we_keep(test_4D01_strings, barcode_4)    
a=1
#endsIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3
#endIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3


#given csv file and desired file name for kept CCSs, filters CCS
def filter_ccs(all_ccs_filename, filtered_filename, barcode):
    retained_csv_list = []
    with open(all_ccs_filename,'r') as f:
        if(os.path.exists(all_ccs_filename)):
            test_var = 1
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                line_list = line.split(',')
                ccs = line_list[1]
                if should_we_keep(ccs, barcode):
                    retained_csv_list.append(line_list)
            try:
                f = open(filtered_filename, 'w+')
                for entry in retained_csv_list:
                    f.write(entry[0] + ',' + entry[1])
                    sleep(0.2)
                    f.write('\n')
                f.close()  
            except:
                print("Exception: something went wrong working with filtered_filename: " + filtered_filename)
                print("here's the traceback: " + traceback + print_exc())
        else:
            print("Sorry!! all_ccs_filepath: " + all_ccs_filename + "doesn't exist!!") 






